from dash import Dash, html

style = {
    'backgroundColor': 'blue',
    'width': '20%'    
    }

def SideBar(children, **kwargs):
    return html.Div(children, style=style)




