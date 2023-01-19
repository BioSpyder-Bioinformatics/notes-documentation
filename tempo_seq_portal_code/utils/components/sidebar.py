import dash
from dash import Dash, html, dcc, Output, Input, State, callback

style = {
    
    }

def SideBar(app, children='This is the sidebar', **kwargs):
    #These need to be styled here as the css can't allow two ids to be the same
    div_large = html.Div(
        [   
            html.Button('Click me, im large', id='sidebarBtn'),
            html.Div([
                dcc.Link(page['name']+" | ", href=page['path']) for page in dash.page_registry.values()
            ])
        ], 
        style=style
    )

    div_small = html.Div(
        [   
            html.Button('small', id='sidebarBtn'),
            html.Div([
                dcc.Link(page['name']+" | ", href=page['path']) for page in dash.page_registry.values()
            ])
        ], 
        style=style
    )
    


    
    #Callback to change the style of navbar
    @app.callback(
        #Output style (width), children to sidebar
        Output(component_id='sidebarDiv', component_property='style'),
        Output(component_id='sidebarDiv', component_property='children'),
        #Output style (width) to main div
        Output(component_id='contentDiv', component_property='style'),
        Input(component_id='sidebarBtn', component_property='n_clicks'),
        State(component_id='sidebarDiv', component_property='style'),
        prevent_initial_call=True
    )
    def change_sidebar_style(click, style):
        print('CLICKED', click, style)
        print(dash.page_registry.values())
        if style['width'] == '12%':
            return {'width': '3%'}, div_small, {'width': '97%'}
        else:
            return {'width': '12%'}, div_large, {'width': '88%'}
        

    return div_large



