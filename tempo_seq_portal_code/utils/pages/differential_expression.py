import dash
from dash import html, dcc

dash.register_page(__name__, path='/differential_expression', name='Differential Expression')

layout = html.Div(
    dcc.Markdown('# DIFFERENTIAL EXPRESSION')
)

