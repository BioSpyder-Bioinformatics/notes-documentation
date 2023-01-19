import dash
from dash import html, dcc

dash.register_page(__name__, path='/sample_sheet', name='Sample Sheet Generator')

layout = html.Div(
    dcc.Markdown('# SAMPLE SHEET')
)





