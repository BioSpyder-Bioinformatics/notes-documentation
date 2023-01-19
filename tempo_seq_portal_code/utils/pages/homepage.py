from dash import Dash, html, dcc
import dash


style = {
    
}


#For later on to map to navbar
dash.register_page(__name__, path='/', name='Homepage')#as it is homepage


#instead of returning things, we just declare the layout
layout = html.Div(
    [dcc.Markdown('# HOMEPAGE'),
    html.Div('Hello there',
    style={
        'height':'300px',
        'width':'90%',
        'backgroundColor':'blue'
    })], 
    style=style,
    className='row')





# def HomePage(children = 'This is the homepage', **kwargs):
#     return html.Div(children, style=style, className='row')

