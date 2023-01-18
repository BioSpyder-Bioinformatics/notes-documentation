import dash
from dash import Dash, html, dcc, Output, Input, State, callback

style = {
    
    }

def SideBar(app, children='This is the sidebar', **kwargs):
    div = html.Div(
        [   
            html.Button('Click me', id='sidebarBtn'),
            html.Div([
                dcc.Link(page['name']+" | ", href=page['path']) for page in dash.page_registry.values()
            ])
        ], 
        style=style
    )
    
    #Callback to change the style of navbar
    @app.callback(
        Output(component_id='sidebarDiv', component_property='style'),
        Input(component_id='sidebarBtn', component_property='n_clicks'),
        State(component_id='sidebarDiv', component_property='style'),
        prevent_initial_call=True
    )
    def change_sidebar_style(click, style):
        print('CLICKED', click, style)
        if not style:
            return {'backgroundColor': 'purple'}
        if style['backgroundColor'] == 'green':
            return {'backgroundColor': 'purple'}
        else:
            return {'backgroundColor': 'green'}
        

    return div



