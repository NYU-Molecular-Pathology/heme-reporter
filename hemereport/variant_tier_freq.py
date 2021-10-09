import pandas as pd
from itertools import chain

def calc_freq(data):
    calc_freq = data.groupby(["variants", "tier","user"]).size().reset_index(name="Times")
    return(calc_freq)

def variant_tier_frequency(input_data):
    data = input_data
    freq_variants = calc_freq(data)
    freq_variants['user'] = "User" + ':' + freq_variants['user'].map(str)
    freq_variants['tier'] = "Tier" + ':' + freq_variants['tier'].map(str)
    freq_variants['Times'] = "Freq" + ':' + freq_variants['Times'].map(str)
    freq_variants_values = freq_variants.groupby(['variants'])[['user','tier','Times']].apply(lambda x: x.values.tolist())
    freq_variants_values_df = freq_variants_values.to_frame()
    freq_variants_values_df.columns = ['Freq']
    freq_variants_values_df["Freq"] = [', '.join(chain.from_iterable(x)) for x in freq_variants_values_df["Freq"]]
    return(freq_variants_values_df)
    






