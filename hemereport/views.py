from . import process_input, variant_tier_freq
from django.shortcuts import render, redirect
from django.urls import reverse
from .models import PMKBdb, AllAbberations, VariantTiering, VariantSpecificComment
from .forms import AllAbberationsForm
import pandas as pd
import numpy as np
import json, math

# Create your views here.
def home(request):
    return render(request, 'hemereport/home.html', {} )

def upload_input_report(request):
    return render(request, 'hemereport/upload_input_report.html', {})

def generate_abberations(request):
    ## importing from csv using pandas https://deeplearning.lipingyang.org/import-csv-using-pandas-to-django-models/ ##
    csv_file = request.FILES["csv_file"]
    csv_file_name = request.FILES['csv_file'].name
    name = csv_file_name.split(".csv")
    filename = name[0]
    if AllAbberations.objects.filter(runid=filename).exists():
        variant_comment = pd.DataFrame(list(VariantSpecificComment.objects.all().values()))
        variant_comment_list = list(variant_comment.variants)
        ### Get samples list and render as dropdown ###
        all_data = pd.DataFrame(list(AllAbberations.objects.filter(runid=filename).values()))
        all_data['DefaultTier'] = all_data['variants'].apply(lambda x: 'Yes' if x in variant_comment_list else 'No')
        ## Add specimen type dropdown list ##
        specimentype= ["peripheral blood","bone marrow aspirate","other"]
        ## Get all the variants and tier from VariantTiering table ##
        tiering_history = pd.DataFrame(list(VariantTiering.objects.all().values()))
        ## if freq table is not empty then calc freq using the exsisting values ##
        if not tiering_history.empty:
        ## Get the variant freq for the selected list from table ##
            freq_variants_values = variant_tier_freq.variant_tier_frequency(tiering_history)
            all_data_freq_merge = pd.merge(all_data,freq_variants_values,how='left',on='variants')
            samples = all_data['sample'].unique()
            #abberations_df = pd.DataFrame(list(AllAbberations.objects.all().values()))
            abberations = all_data_freq_merge.drop_duplicates()
            json_records = abberations.reset_index().to_json(orient ='records') 
            data_df = []
            data_df = json.loads(json_records)
        else:
            samples = all_data['sample'].unique()
            abberations = all_data.drop_duplicates()
            json_records = abberations.reset_index().to_json(orient ='records') 
            data_df = []
            data_df = json.loads(json_records)

        return render(request, 'hemereport/generate_abberations.html', {'samples':samples, 'data_df':data_df,'specimentype':specimentype}) 
    else:
        input_df = pd.read_csv(csv_file)
        input_df['Runid'] = filename
        ## Strip of all samples except TMs ##
        input_df_TMonly = input_df.loc[(input_df['Sample'].str.startswith("TM",na=False))]
        input_df_TMonly_setindex = input_df_TMonly.reset_index(drop=True)
        ## Get the indexes/row numbers where gene is FLT3 and not SNV types ##
        idx = np.where((input_df_TMonly_setindex['Genes']=="FLT3") & (input_df_TMonly_setindex['Type'] != "SNV"))
        ## Convert numpy array to list of lists ##
        lst = np.array(idx).tolist()
        ## Convert list of lists to flatlist ##
        flat_list = [item for sublist in lst for item in sublist]
        ## for each of those rows with FLT3 INDELS/ITDs
        for index in flat_list:
            input_df_TMonly_setindex.loc[index,'Type'] = input_df_TMonly_setindex.loc[index,'Type']+" "+input_df_TMonly_setindex.loc[index,'Length'].astype(str).strip(".0")+"(bp)"
        input_df_TMonly_setindex['Genes'] = input_df_TMonly_setindex['Genes'].fillna("NEGATIVE")
        input_df_cols = input_df_TMonly_setindex[['Runid','Sample','Locus','Genes','Type','Exon','Transcript','Coding','Variant.Effect','Length','Frequency','AA','Coverage']]

        myeloid_report = input_df_cols
        PMKB = pd.DataFrame(list(PMKBdb.objects.all().values()))

        input_all_abberations = process_input.main(myeloid_report,PMKB)
        variant_comment = pd.DataFrame(list(VariantSpecificComment.objects.all().values()))
        variant_comment_list = list(variant_comment.variants)

        ## If there are any variant specific comment genes present update their tier, interpretations and citations in all abberations table  ##
        for variant in variant_comment.variants:
            if variant in input_all_abberations.values:
                rows_found = input_all_abberations[input_all_abberations.Variants == variant].index
                variants_comment = variant_comment[variant_comment.variants == variant]
                variants_tier_toupdate = variants_comment['tier']
                variants_interpret_toupdate = variants_comment['interpretations']
                variants_cit_toupdate = variants_comment['citations']
                for indexes in rows_found:
                    input_all_abberations.at[indexes,'tier'] = variants_tier_toupdate.values[0]
                    input_all_abberations.at[indexes,'interpretations'] = variants_interpret_toupdate.values[0]
                    input_all_abberations.at[indexes,'citations'] = variants_cit_toupdate.values[0]
        ## Else use the original/default all abberations table ##
            else:
                input_all_abberations = input_all_abberations

        input_all_abberations_final = input_all_abberations.drop_duplicates()
        row_iter_all = input_all_abberations_final.iterrows()

        objs_all = [
            AllAbberations(
                    runid = row['Runid'],
                    sample = row['Sample'],
                    genes = row['genes_myeloid'],
                    variant_type = row['Type'],
                    transcript = row['Transcript'],
                    variants = row['Variants'],
                    length_bp = row['Length'],
                    vaf = row['Frequenc'],
                    exon = row['Exon'],
                    amino_acid_change = row['AA'],
                    coverage = row['Coverage'],
                    tier = row['tier'],
                    locus = row['Locus'],
                    interpretations = row['interpretations'],
                    citations = row['citations'],
                    variant_pmkb = row['variant'],
                )
                for index, row in row_iter_all
                ]

        AllAbberations.objects.bulk_create(objs_all)

        ### Get samples list and render as dropdown ###
        all_data = pd.DataFrame(list(AllAbberations.objects.filter(runid=filename).values()))
        all_data['DefaultTier'] = all_data['variants'].apply(lambda x: 'Yes' if x in variant_comment_list else 'No')

        #all_data_final = all_data[all_data['variant_pmkb'].notna()]
        specimentype= ["peripheral blood","bone marrow aspirate","other"]
        ## Get all the variants and tier from VariantTiering table ##
        tiering_history = pd.DataFrame(list(VariantTiering.objects.all().values()))
        if not tiering_history.empty:
            ## Get the variant freq for the selected list from table ##
            freq_variants_values = variant_tier_freq.variant_tier_frequency(tiering_history)
            all_data_freq_merge = pd.merge(all_data,freq_variants_values,how='left',on='variants')
            samples = all_data['sample'].unique()
            abberations = all_data_freq_merge.drop_duplicates()
            json_records = abberations.reset_index().to_json(orient ='records') 
            data_df = []
            data_df = json.loads(json_records)
        else:
            samples = all_data['sample'].unique()
            abberations = all_data.drop_duplicates()
            json_records = abberations.reset_index().to_json(orient ='records') 
            data_df = []
            data_df = json.loads(json_records)
        return render(request, 'hemereport/generate_abberations.html', {'samples':samples, 'data_df':data_df,'specimentype':specimentype})

## can go by runid and sample, but leads to that sample info only ##
#def index(request, runid, sample):
def index(request, runid): 
    #allabberations_df = pd.DataFrame(list(AllAbberations.objects.all().values()))
    specimentype= ["peripheral blood","bone marrow aspirate","other"]
    allabberations_df = pd.DataFrame(list(AllAbberations.objects.filter(runid=runid).values()))
    tiering_history = pd.DataFrame(list(VariantTiering.objects.all().values()))
    if not tiering_history.empty:
        ## Get the variant freq for the selected list from table ##
        freq_variants_values = variant_tier_freq.variant_tier_frequency(tiering_history)
        allabberations_merge = pd.merge(allabberations_df,freq_variants_values,how='left',on='variants')
        json_records = allabberations_merge.reset_index().to_json(orient ='records') 
        allabberations = []
        allabberations = json.loads(json_records)
        samples = allabberations_df['sample'].unique()
    else:
        allabberations = AllAbberations.objects.all()
        samples = allabberations_df['sample'].unique()
    return render(request,"hemereport/show.html",{'allabberations':allabberations,'samples':samples,'specimentype':specimentype}) 

def index_runid_sample(request, runid, sample):
    specimentype= ["peripheral blood","bone marrow aspirate","other"]
    allabberations_df = pd.DataFrame(list(AllAbberations.objects.filter(runid=runid,sample=sample).values()))
    tiering_history = pd.DataFrame(list(VariantTiering.objects.all().values()))
    if not tiering_history.empty:
        ## Get the variant freq for the selected list from table ##
        freq_variants_values = variant_tier_freq.variant_tier_frequency(tiering_history)
        allabberations_merge = pd.merge(allabberations_df,freq_variants_values,how='left',on='variants')
        json_records = allabberations_merge.reset_index().to_json(orient ='records') 
        allabberations = []
        allabberations = json.loads(json_records)
        selected_sample = sample
        runid = runid
    else:
        allabberations = AllAbberations.objects.filter(runid=runid,sample=sample)
        selected_sample = sample
        runid = runid
    return render(request,"hemereport/show_persample.html",{'allabberations':allabberations,'sample':selected_sample, 'runid':runid,'specimentype':specimentype}) 

def edit(request, id):
    allabberations = AllAbberations.objects.get(id=id)
    return render(request, 'hemereport/edit.html', {'allabberations': allabberations})

def update(request, id):
    allabberations_id = AllAbberations.objects.get(id=id)
    runid = allabberations_id.runid
    sample = allabberations_id.sample
    allabberations = AllAbberations.objects.filter(id=allabberations_id.id,runid=runid).first()
    form = AllAbberationsForm(request.POST, instance = allabberations)
    if form.is_valid():
        form.save()
        #return redirect('index',runid=runid)
        return redirect('index_runid_sample',runid=runid,sample=sample)
    return render(request, 'hemereport/edit.html', {'allabberations': allabberations} ) 

def preview_report(request):
    selected_variants = request.POST.getlist('selected_values')
    specimentype = request.POST['specimentype']
    if not selected_variants:
        sample = request.POST['Sample']
        if "Q" in sample:
            tumor_sample = sample.split("-Q")[0]
        else: 
            tumor_sample = sample.split("-B")[0]
        ## Get the run id ##
        df = pd.DataFrame(list(AllAbberations.objects.filter(sample=sample).values()))
        ## Get the runid to redirect to index page with runid ##
        variants_runid = df.runid.unique()
        runid = variants_runid[0]
        return render(request, 'hemereport/negative_report.html', {'tumor_sample': tumor_sample, 'runid': runid, 'specimentype': specimentype})
    else:
        sample = request.POST['Sample']
        if "Q" in sample:
            tumor_sample = sample.split("-Q")[0]
        else: 
            tumor_sample = sample.split("-B")[0]
        ## Get selected checkbox variants ##
        selected_variants = request.POST.getlist('selected_values')
        res = [sub.replace('/', '') for sub in selected_variants]
        ## query the abberations if id is present is res list ##
        selected_variants_df = pd.DataFrame(list(AllAbberations.objects.filter(id__in=res).values()))
        ## Trying to save selected variants into tiering table for freq calculation per user##
        selected_variants_df['user'] = request.user.username
        selected_variants_df['sample'] = tumor_sample
        ## Get the runid to redirect to index page with runid ##
        variants_runid = selected_variants_df.runid.unique()
        runid = variants_runid[0]

        tiers_to_update = selected_variants_df['tier']
        user_to_update = request.user.username

        ## If sample already exsists in the variant tiering table, then update tier for that user and id ##
        get_variant_id = VariantTiering.objects.filter(abberation_id__in=res,runid=runid,sample=tumor_sample).values('id')
        print(get_variant_id)
        ## if it exsists then update the tier for that id ##
        if get_variant_id.exists():
            tier_user_updated = VariantTiering.objects.filter(id__in=get_variant_id)
            for variant_tier_id, tiers_to_update_value in zip(tier_user_updated,tiers_to_update):
                variant_tier_id.tier = tiers_to_update_value
                variant_tier_id.user = user_to_update
                variant_tier_id.save()
        ## Else create the whole set of objects and import into variant tiering table ##
        else:
            row_iter_selected = selected_variants_df.iterrows()

            objs_all_selectedvariants = [
            
            VariantTiering(
                    abberation_id = row['id'],
                    user = row['user'],
                    runid = row['runid'],
                    sample = row['sample'],
                    genes = row['genes'],
                    variants = row['variants'],
                    tier = row['tier'],
                )
                for index, row in row_iter_selected
            ]
            VariantTiering.objects.bulk_create(objs_all_selectedvariants)
        
        
        ## Need to categorize df based on DNA results,RNA results and Tiers ##
        variants_table = selected_variants_df[['genes','variants','tier','variant_type','vaf','coverage','transcript','locus','exon','length_bp','interpretations','citations']]
        ## DNA results ##
        DNA_results = variants_table[(variants_table['variant_type'] != 'FUSION')]
        DNA_results_selected_cols = DNA_results[['genes','variants','tier','variant_type','vaf','coverage','transcript','locus','exon']]
        DNA_results_selected_cols_final = DNA_results_selected_cols.copy()
        ## Need only coding and protein change in table for Variant ##
        DNA_results_selected_cols_final['variants'] = DNA_results_selected_cols_final['variants'].str.split(", ").str[1:3].str.join(', ')
        DNA_results_selected_cols_df = DNA_results_selected_cols_final.sort_values('tier')
        DNA_results_final = DNA_results_selected_cols_df.to_dict('records')

        ## DNA Tiers ##
        DNA_Tier1 = DNA_results[(DNA_results['tier'] == 1)]
        DNA_Tier1_final= DNA_Tier1.to_dict('records')
        DNA_Tier2 = DNA_results[(DNA_results['tier'] == 2)]
        DNA_Tier2_final= DNA_Tier2.to_dict('records')
        DNA_Tier3 = DNA_results[(DNA_results['tier'] == 3)]
        DNA_Tier3_final = DNA_Tier3.to_dict('records')

        ## If variant specific comments are present then add it to render (currently only KIT in tier3 provided by TG ) ##
        variant_comment = pd.DataFrame(list(VariantSpecificComment.objects.all().values()))

        variant_comment_tier3 = pd.DataFrame()
        for variant in variant_comment.variants:
            if variant in DNA_Tier3.values:
                variant_comment_tier3 = variant_comment_tier3.append(variant_comment)
            else:
                variant_comment_tier3 = variant_comment_tier3.append(pd.Series(), ignore_index=True)
        
        ## RNA results ##
        RNA_results = variants_table[(variants_table['variant_type'] == 'FUSION')]
        RNA_results_final = RNA_results.to_dict('records')

        all_values = {'tumor_sample': tumor_sample,'DNA_results_final': DNA_results_final, 'DNA_Tier1_final':DNA_Tier1_final, 'DNA_Tier2_final':DNA_Tier2_final, 'DNA_Tier3_final': DNA_Tier3_final, 'variant_comment_tier3':variant_comment_tier3, 'RNA_results_final': RNA_results_final, 'runid': runid, 'specimentype': specimentype}
        return render(request, 'hemereport/insert_selected_variants.html', all_values)
