import pandas as pd
from pathlib import Path
from pandas import DataFrame

match = ["ENERGY", "dG"] #find header rows in ct file
def convert_ct_to_dataframe(file):
    """
        Convert the ct file to a dataframe. This also gives the number of structures in the ct file.
        This method should be inside of a with expression.
        :param filename: the file to read and convert into a dataframe.
        :return: (Dataframe, int structure_count)
        """

    # dtype={"baseno": int, "base": str, "bs_bind": int}
    file.seek(0)
    ct_df =  pd.read_csv(file, sep='\\s+', usecols=[0,1,4], names=["baseno", "base", "bs_bind"], engine='python')
    initial_row_count = len(ct_df)
    ct_df = ct_df[~ct_df["base"].isin(match)]
    structure_count = initial_row_count - len(ct_df)
    ct_df = ct_df.astype({"baseno": int, "base": str, "bs_bind": int})
    return (ct_df, structure_count)

def getSSCountDF(ct_dataframe : DataFrame, save_to_file: bool = None, output_file: Path = None) -> DataFrame:
    """
    Get the SSCount from a dataframe
    :param ct_dataframe: ct converted to a dataframe
    :return: Dataframe
    """

    df_grouped = ct_dataframe.groupby(['baseno', 'base'], as_index = False).agg(lambda x: x[x == 0].count())
    df_grouped.rename(columns={'bs_bind': 'sscount'}, inplace=True)
    df_grouped = df_grouped.reindex(columns=['baseno','sscount','base'])


    df_grouped['base'] = df_grouped.base.replace('T', 'U')

    if save_to_file: df_grouped.to_csv(output_file, index=False, header=False)

    return df_grouped