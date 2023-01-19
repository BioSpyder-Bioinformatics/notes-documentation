import dash
from dash import html, dcc

dash.register_page(__name__, path='/aligner', name='Aligner')

layout = html.Div(
    dcc.Markdown('# ALIGNER')
)

