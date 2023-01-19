import dash
from dash import html, dcc

dash.register_page(__name__, path='/support', name='Support')

layout = html.Div(
    dcc.Markdown('# SUPPORT PAGE')
)

