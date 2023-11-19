import pandas as pd
import numpy as np


def df_col_to_dict(df, values_column:str, names_column=None):
    '''
    create a dict for rownames as keys and column is values.
    if names_column is None, will use the index as keys.
    '''
    if names_column is None:
        return dict(zip(df.index, df[values_column]))
    return dict(zip(df[names_column], df[values_column]))
    
    
def merge_col(df1, df2, values_column, by.x=None, by.y=None, suffixes=['x','y'],inplace=False):
    '''
    values_column - name of column from df2 that should be added to df1.    
    by.x, by.y - if None - will use the index, otherwise should be a column name (str).
    '''
    merged = df1 if inplace else df1.copy()
    if values_column not in df1.columns:
        name = values_column
    else:
        name = values_column + '_' + suffixes[1]   
        merged.rename(columns={values_column:values_column + '_' + suffixes[0]}, inplace=True)
    translation = df_col_to_dict(df2, values_column, names_column=by.y)    
    if by.x is None:
        return merged[name] = merged.index.map(translation)
    return merged[name] = merged[by.x].map(translation)


def pdist1(df, column1, column2, type="normal", int_type='float32', diagonal_na=True):
    '''type - "normal" for cell by cell, "long" for three columns (cell1, cell2, distance)'''
    from scipy.spatial.distance import pdist, squareform
    square_distance_matrix = squareform(pdist(df[[column1, column2]].values, metric='euclidean') ) # cell by cell distances matrix
    square_distance_matrix = np.round(square_distance_matrix).astype(int_type)
    if diagonal_na:
        np.fill_diagonal(square_distance_matrix, np.nan)
    distance_df = pd.DataFrame(square_distance_matrix, columns=df.index, index=df.index)
    if type=='long':
        distance_df = distance_df.where(np.triu(np.ones(distance_df.shape)).astype(bool))
        distance_df = pivot_longer_matrix(distance_df)
    return distance_df
    
    
def pivot_longer_matrix(df, colnames=['row','column','value']):
    '''transforms a cell by cell dataframe into three columns (row, column, value)'''
    df = df.stack()
    df = df.rename_axis(["row","column"]).reset_index()
    df.columns = colnames
    return df
    
def memory(df, every_col=False):
    if every_col:
        for col in df.columns:
            print(f'{col}\t\t{df[col].memory_usage()/1000000000} GB')
        print('- - - - - - - - - -')
    print('total:\t\t' + str(df.memory_usage().sum()/1000000000) + ' GB')
    
    
    
    
    
    




# general pandas functions
pd.set_option('display.max_columns', 99)

df = pd.read_csv(path, index_col='cell_id')
df = pd.DataFrame({'col1':[1,2], 'col2':[3,4]} # df from dictionary
df.shape #(rows, columns)

df['col1'].value_counts() # summary (count table of a column)
df.info() 


# types
column_type = df['your_column'].dtype # check type
df['your_column'] = df['your_column'].astype('int32') # change type

#rename columns
df.columns = ['a', 'b', 'c']
df.columns = [x.lower() for x in df.columns] # all to lower
df.columns = df.columns.str.replace('-','_') # replace '-' with '_' (grepl)
df = df.rename(columns={'x': 'a', 'y': 'b'})

# index operations
df = df.set_index('cell_id') # set index
df = df.reset_index()
df = df.sort_index() # ascending=True

#subset columns
your_column = df['your_column']
df_subset = df[['your_column1', 'your_column2']]
df.loc[:, df.columns.str.contains('x')]  # grepl colnames()
df.loc[df['col1'].str.contains('x', na=False),['col1','col2']]

#subset rows
df_subset = df.loc[df['your_column'] > 10] # based on condition
df_subset = df.loc[df['your_column'].isin(['val1','val2','val3'])]
df_subset = df.loc[['index_name1', 'index_name2']] # based on index name
df_subset = df.iloc[0:5] # based on  row number
df_subset = df.iloc[[4,7,42]] # based on row number

# subset both rows and columns
df.iloc[[0,5],4] # rows, then columns
df.iloc[0:5,[1,2,3]]
df.loc[['row1','row42'],'col1':'col42']
df.iloc[0:5, df.columns.get_indexer(['col1', 'col42'])]
df_subset = df.loc[df['col1'] > 10, 'col1':'col42']
df_subset = df.loc[df['col1'] > 10 & df['col1'] < 0, 'col1':'col42'] # |,&
df_subset = df.loc[~df['col1'] > 10 ,:] # |,&

df.at[row_label, column_label] # faster than .loc for single value

# row operations 
df['col42'].apply(foo)
df['col42'].apply(lambda x: np.log2(x))
df['col42'].apply(lambda x: x.lower())

df.apply(sum) # sum of all columns. apply the function on all rows for each column
df.apply(sum, axis='columns') # sum of all rows. apply the function on all columns for each row
df.apply(pd.Series.min) # minimal value of each row

df.applymap(min) # minimal value from the whole dataframe
df.applymap(np.log10) # log10 on the whole dataframe

df['col42'].map({'val1':'sub1', 'val2':'sub2'}) # change values based on dictionary. all other will be NaN
df['col42'].replace({'val1':'sub1', 'val2':'sub2'}) # change values based on dictionary. all other will be left untouched


# add/remove
df['combined'] = df['col1'] + '_' + df['col2'] # merge two columns
df[['col1','col2']] = df['col42'].str.split(' ', expand=True) # split a column
df = df.drop(columns=['col1','col2']) # delete columns
df = df.drop(index=['row2','col1']) # delete columns

df.append({'col1': 'val42'}, ignore_index=True)

# append rows
df = df.append(df2, ignore_index=True)
# delete rows
filt = df['col42'] == 'X'
df.drop(index=df[filt].index)

# merge
merge(df1, df2, by="gene", all=TRUE)                -   pd.merge(df1, df2, on="gene", how="outer")
merge(df1, df2, by="gene", all=FALSE)               -   pd.merge(df1, df2, on="gene", how="inner")
merge(df1, df2, by="gene", all.x=TRUE, all.y=FALSE) -   pd.merge(df1, df2, on="gene", how="left")      
merge(df1, df2, by="gene", all.x=FALSE, all.y=TRUE) -   pd.merge(df1, df2, on="gene", how="right")   
