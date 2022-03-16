
def get_number_and_category_cols(parsed_df, y_label):

    category_columns = list(parsed_df.select_dtypes(include='object').columns)
    category_columns.remove("ID")
    data_df = parsed_df
    data_df.drop("ID", axis=1)

    for col in category_columns:
        data_df[col] = data_df[col].astype(str)

    number_columns = list(parsed_df.select_dtypes(include='float').columns)
    number_columns.remove(y_label)

    return number_columns, category_columns
