import pandas as pd

# Create a sample DataFrame
data = {'Column_A': [1, 2, 3, 4],
        'Column_B': [5, 6, 7, 8],
        'Column_C': [9, 10, 11, 12],
        'Column_D': [13, 14, 15, 16]}

df = pd.DataFrame(data)

# Get the index position of the second column (index starts from 0)
column_to_multiply = 1

# Perform elementwise multiplication with broadcasting
result = df.iloc[:, [column_to_multiply]].values * df.iloc[:, 3:].values

# Create a new DataFrame with the result
result_df = pd.DataFrame(result, columns=df.columns[3:])

# Display the result
print(result_df)