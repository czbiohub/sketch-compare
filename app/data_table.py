import dash
import dash_table

def create_datatable(df):
    dash_table.DataTable(
    #     id='table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
    )
    
    return dash_table.DataTable
