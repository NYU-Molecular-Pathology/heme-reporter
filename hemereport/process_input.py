import pandas as pd
import re

def process_files(myeloid_report_df,PMKB,selected_sample,type):
    myeloid_values = pd.DataFrame()
    pmkb_values = pd.DataFrame()

    ### For each sample get the required columns from input report ###
    sample_filter_df = myeloid_report_df[(myeloid_report_df['Sample'] == selected_sample)]
    ## split Genes column with multiple gene names as separate rows (https://riptutorial.com/pandas/example/25462/split--reshape--csv-strings-in-columns-into-multiple-rows--having-one-element-per-row) ##
    ## For DNA samples example: TET2,TET2-AS1 to TET2 and TET2-AS1 as separate rows ##
    reshaped_df = sample_filter_df.set_index(sample_filter_df.columns.drop('Genes',1).tolist()).Genes.str.split(',', expand=True).stack().reset_index().rename(columns={0:'Genes'}).loc[:,sample_filter_df.columns]
    
    reshaped_df['Gene'] = reshaped_df['Genes'] ## to get fusion gene original name

    ## For RNA samples example: BCR(2)-ABL1(4) to BCR(2) and ABL1(4) as separate rows ##
    reshaped_df_final = reshaped_df.set_index(reshaped_df.columns.drop('Genes',1).tolist()).Genes.str.split(" - ", expand=True).stack().reset_index().rename(columns={0:'Genes'}).loc[:,reshaped_df.columns]

    ## Add Variants column to get myeloid report values ##
    if type == "DNA":
        reshaped_df['Variants'] = reshaped_df[['Genes','Coding','AA']].apply(lambda x: ', '.join(x), axis=1)
        myeloid_values = myeloid_values.append(reshaped_df)
    else:
        reshaped_df_RNA = reshaped_df_final
        reshaped_df_RNA['Variants'] = "NA" # No need to combine variants (no AA change)
        reshaped_df_RNA['Genes'] = reshaped_df_RNA['Genes'].str.replace(r" ?\([^)]+\)","",regex=True) ## Remove exon numbers and brackets from gene names
        myeloid_values = myeloid_values.append(reshaped_df_RNA)

    for gene in reshaped_df_final.Genes.unique():
        if type == "DNA":
            myeloid_gene_in_pmkb = PMKB[(PMKB['gene'] == gene)]
            myeloid_gene_notin_pmkb = reshaped_df_final[(reshaped_df_final['Genes'] == gene)]
            merged = myeloid_gene_notin_pmkb.merge(myeloid_gene_in_pmkb,left_on='Gene',right_on='gene',how='left')
            pmkb_values = pmkb_values.append(merged)
        else:
            remove_exons_gene = re.sub(r" ?\([^)]+\)", "", gene)
            gene_final = remove_exons_gene
            myeloid_gene_in_pmkb = PMKB[(PMKB['gene'].str.contains(gene_final)) & (PMKB['variant'].str.contains('fusion'))]
            myeloid_gene_notin_pmkb = reshaped_df_RNA[(reshaped_df_RNA['Genes'] == gene_final) & (reshaped_df_RNA['Type'] == "FUSION")]
            merged = pd.merge(myeloid_gene_notin_pmkb,myeloid_gene_in_pmkb,how='left',left_on='Genes',right_on='gene')
            pmkb_values= pmkb_values.append(merged)

    #pmkb_values_nona = pmkb_values[pmkb_values['variant'].notna()]
    pmkb_valuesfinal = pmkb_values.drop_duplicates(subset=['Gene','Type','interpretations','Locus'],keep="last")
    myeloid_valuesfinal = myeloid_values.drop_duplicates(subset=['Genes','Type','Variants','Locus'], keep="last")

    merged_table_dna = pd.merge(myeloid_valuesfinal,pmkb_valuesfinal,how='left',on='Genes')
    merged_table_rna = pd.merge(myeloid_valuesfinal,pmkb_valuesfinal,how='left',on='Genes')
    return merged_table_dna, merged_table_rna   
        
def all_results(myeloid_report,PMKB,selected_sample):
    myeloid_report_persample = myeloid_report[(myeloid_report['Sample'] == selected_sample)]
    ## If Sample has SNV,INDEL then combine DNA+RNA 
    myeloid_report_persample_type = myeloid_report_persample.Type.unique()
    type_status = any(x in myeloid_report_persample_type for x in ['SNV', 'INDEL'])
    ## For negative cases, if Genes is "NEGATIVE" then just return Run and sample info ##
    ## Checking if Genes contians NEGATIVE ##
    negative_found = myeloid_report_persample[myeloid_report_persample['Genes'].str.contains('NEGATIVE')]
    if type_status == True: 
        ## Get both DNA and RNA results in other cases which have both info ##
        ### Get DNA results table ###
        myeloid_report_df_DNA = myeloid_report[(myeloid_report['Type'] != 'FUSION')  & (myeloid_report['Sample'] == selected_sample)]
        ## need to include FLT3 indels and have only FLT3ITDs ##
        ### Get unique sample list from myeloid report ###
        DNA_results_processed = process_files(myeloid_report_df_DNA,PMKB,selected_sample,"DNA")[0]
        DNA_results_final = DNA_results_processed[['Runid_x','Sample_x','Genes','Type_x','Transcript_x','Variants','Length_x','Frequency_x','Exon_x','AA_x','Coverage_x','tier','Locus_x','interpretations','citations','variant']]
        DNA_results_final.columns = DNA_results_final.columns.str.rstrip('_x|_y')
        DNA_results_final.rename(columns={'Genes':'genes_myeloid'}, inplace=True)
        DNA_results_final_tier = DNA_results_final.copy()
        DNA_results_final_tier['tier'] = DNA_results_final_tier['tier'].fillna(0)
        DNA_results_final_df = DNA_results_final_tier.sort_values('tier')
        DNA_results = DNA_results_final_df.drop_duplicates()

        myeloid_report_df_RNA = myeloid_report[(myeloid_report['Type'] == 'FUSION') & (myeloid_report['Sample'] == selected_sample)]
        if myeloid_report_df_RNA.empty:
            DNA_RNA = DNA_results
            abberations_all = DNA_RNA.drop_duplicates()
            return abberations_all
        else:
            RNA_results_processed = process_files(myeloid_report,PMKB,selected_sample,"RNA")[1]
            RNA_results_final = RNA_results_processed[['Runid_x','Sample_x','Gene_x','Type_x','Transcript_x','Variants_y','Length_x','Frequency_x','Exon_x','AA_x','Coverage_x','tier','Locus_x','interpretations','citations','variant']]
            RNA_results_final.columns = RNA_results_final.columns.str.rstrip('_x|_y')
            RNA_results_final.rename(columns={'Gene':'genes_myeloid'}, inplace=True)
            RNA_results_final_df = RNA_results_final.sort_values('tier')
            RNA_results = RNA_results_final_df.drop_duplicates()
            RNA_results_nona = RNA_results[RNA_results['genes_myeloid'].notna()]

            DNA_RNA = DNA_results.append(RNA_results_nona)
            abberations_all = DNA_RNA.drop_duplicates()
            return abberations_all
        
    ## If negative cases are found return minimal info ##
    elif negative_found['Genes'].count() > 0:
        runid = myeloid_report_persample.Runid.unique()
        myeloid_report_persample_final = pd.DataFrame(columns=['Runid','Sample','Genes','Type','Transcript','Variants','Length','Frequency','Exon','AA','Coverage','tier','Locus','interpretations','citations','variant'])
        myeloid_report_persample_final = myeloid_report_persample_final.append({'Runid': runid[0], 'Sample': selected_sample, 'Genes': "negative",'Type' : "negative",'Transcript': "negative",'Variants': "negative",'Frequency': 0.0,'Exon':0,'AA': "negative",'Coverage': 0,'tier': 0,'Locus' : "negative",'interpretations': "negative",'citations': "negative",'variant': "negative"}, ignore_index=True)
        myeloid_report_persample_final.rename(columns={'Genes':'genes_myeloid'}, inplace=True)
        return myeloid_report_persample_final
    
    ## Else return only FUSION information ##
    else:
        ## If it is only FUSION for a sample, process for RNA results
        RNA_results_processed = process_files(myeloid_report,PMKB,selected_sample,"RNA")[1]
        RNA_results_final = RNA_results_processed[['Runid_x','Sample_x','Gene_x','Type_x','Transcript_x','Variants_y','Length_x','Frequency_x','Exon_x','AA_x','Coverage_x','tier','Locus_x','interpretations','citations','variant']]
        RNA_results_final.columns = RNA_results_final.columns.str.rstrip('_x|_y')
        RNA_results_final.rename(columns={'Gene':'genes_myeloid'}, inplace=True)
        RNA_results_final_df = RNA_results_final.sort_values('tier')
        RNA_results = RNA_results_final_df.drop_duplicates()
        return RNA_results

def main(myeloid_report,PMKB):
    all_abberations = pd.DataFrame()
    samples = myeloid_report['Sample'].unique()
    for sample in samples:
        abberation_per_sample = all_results(myeloid_report,PMKB,sample)
        all_abberations = all_abberations.append(abberation_per_sample)           
    input_all_abberations = all_abberations[['Runid','Sample','genes_myeloid','Type','Transcript','Variants','Length','Frequenc','Exon','AA','Coverage','tier','Locus','interpretations','citations','variant']]
    input_all_abberations_dropna = input_all_abberations.dropna(subset=('Variants','interpretations','citations','variant'), how='all')
    input_all_abberations_tiering = input_all_abberations_dropna.copy()
    input_all_abberations_tiering['tier'] = input_all_abberations_tiering['tier'].fillna(0)
    input_all_abberations_final = input_all_abberations_tiering.reset_index(drop=True)
    input_all_abberations_final = input_all_abberations_final[input_all_abberations_final['genes_myeloid'].notna()]
    input_all_abberations_final_NAreplaced = input_all_abberations_final.fillna(0)
    return input_all_abberations_final_NAreplaced

