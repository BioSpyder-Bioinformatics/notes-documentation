from dash import Dash, html

style = {
    'backgroundColor': 'red',
    'height': '100%',
    'width': '80%'
    
}


def HomePage(children, **kwargs):
    return html.Div(children, style=style, className='row')

