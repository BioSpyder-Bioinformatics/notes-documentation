import dash
from dash import html, dcc

dash.register_page(__name__, path='/qc_tools', name='Quality Control')

layout = html.Div(
    dcc.Markdown('# QUALITY CONTROL')
)

