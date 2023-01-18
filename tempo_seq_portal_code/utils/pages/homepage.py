from dash import Dash, html, dcc
import dash


style = {
    'backgroundColor': 'red',
    'height': '100%',
    'width': '80%'
    
}


#For later on to map to navbar
dash.register_page(__name__, path='/')#as it is homepage


#instead of returning things, we just declare the layout
layout = html.Div(
    dcc.Markdown('#HOMEPAGE'), 
    style=style,
    className='row')





# def HomePage(children = 'This is the homepage', **kwargs):
#     return html.Div(children, style=style, className='row')

