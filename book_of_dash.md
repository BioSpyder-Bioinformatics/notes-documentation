# Python keyword revision
- Class: A blueprint to create object. The class defines attributes (data) and methods (functions) of the objects
- Object: A piece of encapsulated data with associated functionality which is built according to a class definition. Also defined as instance of a class.
- Instantiation: Process of creating an object of a class
- Method: Function associated with object
- Attribute: Variable associated with class or instance
- Class attribute: Variable created staticaly in class definition (aka static attribute)
- Dynamic attribute: Object attribute defined dynamically during program execution
- Instance attribute: Variable that holds data belonging ot only a single object
- Inheritance: Concept that allows you to create new classes as modifications of existing classes by reusing some or all the data in new class



# Important points
- In general we use Plotly Express for creating graphs, but we can also use Plotly Graph Objects which is a lower level interface that is more complicated but allows for better personalisation
- There are two types of components
	+ dash.html.components (HTML) -> structural elements such as headings and dividers to style and position elements on the page
	+ dash.core.components (DCC) -> for core app functionality, includes input fields and figures
	+ Actually there are the bootstrap components too but I'd rather not rely on external libraries for styling
- Variables and DFs declared outside a function are global! It's better to store data in dash core component Storage!!!
- The layout of the app refers to the alignment of the components. The style refers to how the element look, also known as *props*
- Dash can render _cytoscape_, a js library for rendering complex graphs. It can also render VTK for 3D graphs
- Depending on the type of data required we can use `pandas_datareader` to download from API and have it in DF format directly

--------------------------------------------------


# App CSS styling
*You can change style* with a stylesheet as shown below or with the style attribute within the components eg `html.Div([...], style={'color':'red', 'backgroundColor':'yellow'}) `
_In Dash_ the dictionary keys for the style are _camelCased_ (instead of hyphened as CSS is)


### External CSS stylesheet
Made by the creator, customise it for your needs: 
https://codepen.io/chriddyp/pen/bWLwgP.css

It describes width and height of the columns and rows on the page

Usage on app
```py
stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=stylesheets)
```


*Useful classes for stylesheet above (honestly looks like grid)*
`className='class_name'`
- row
	+ If assigned to div, puts all elements in the div on the same row
- n columns
	+ For n between one and twelve (written in letters). Allows to decide how many 'columns' should this component occupy out of 12 (do not exceed 12 columns total as each occupies 1/12% of monitor) `className = 'four columns'`





--------------------------------------------------



# Useful components
Many more in https://dash.plotly.com/dash-core-components


## Useful HTML components
*For all HTML components* that accept children props: the first positional argument is children (taking an array of objects) (`children=[]`); Children props are displayed on the page, if it's text or other components (ex. `html.H1('Hello World!')`-> will display `<h1>Hello world</h1>`)

- html.Div()
	+ Essential component for containing other components. By default will take up the full witdth of the parent container. Consequently the first div you declare takes on the full width of the page. You can leverage this with a flexbox for layout
- html.H1-6()
	+ Headers (1 to 6)
- html.A()
	+ Anchor link for hyperlinked elements
		* 4 props `html.A(id='myId', children='Click here', href='https://google.com', target='_blank'`
			- id -> id for CSS reference
			- children -> regular children prop
			- href -> link destination
			- target -> if `_self` the link opens in the same tab, if `_blank` in a new one




## Useful DCC (Dash Core Components)
Prebuilt components to interact with the app

- dcc.Graph(`id='x', figure={}`)
	+ Component that allows to incorporate data visualisation from Plotly
	+ Main props are `id` and `figure`
		* id -> as always for CSS reference
		* figure -> placeholder for the Plotly chart (placeholder object is empty dict)


- dcc.Dropdown(`id='x', multi=True, options=[], value=[]`)
	+ Component that allows users to choose options from a dropdown menu to dynamically filter data and update graphs
	+ Main props are id, multi, options and value
		* id
		* multi -> allows to choose whether the user can select multiple values at once
		* options -> represents the values the user can choose from when they click the dropdown. It takes in a list of dictionaries in the format {'label': 'Displayed value', 'value': 'reference value'}. Ex with list comprehension: options=[{'label':x, 'value':y} for x,y in zip(labels, values)]
		* value -> represents the default value the dropdown will take on boot


- dcc.RangeSlider(`id='x', min=2000, max=2010, step=1, value=[2002, 2004], marks={2005:'2005', 2006:'06', 2007:'07', ...}`)
	+ Component that allows to select a range of values instead of discrete only
	+ Main props are id, min, max, step, value, marks
		* id
		* min/max -> minimum and maximum value displayed on the slider
		* step -> unit of increment
		* value -> list of currently selected values range
		* marks -> label per value in dict form
			- Extra
				+ allowCross=True (default:false) -> Allows slider handles to cross each other. Set true just for user experience, doesn't change overall functionality


*!!!!!!!!!!!!!*
- dcc.Store(`id='x', storage_type='local', data={}`)
	+ Component used to save dashboard data in memory, so that it can be recalled quickly and efficiently between callbacks. Store is invisible but has to still be declared in the layout section. Limited to 2MB in mobile and 5 to 10MB in desktop only application.
	+ Main props are id, storage_type and data
		* id
		* storage_type -> how we want to store the data: 'session', 'local' or 'memory'. 'session' retains data until the browser tab is closed; 'local' saves the data to the broser until all browing history and cookies are deleted; 'memory' reset data when browser is refreshed.
		* data -> can be dict, list, integer, string or bool. It is not necessary to declare it on instantiation.


- dcc.Interval(`id='x', interval=100, n_intervals=0`)
	+ Component used to automatically update the app without having to refresh the page manually. Typically used with apps that use real-time data
	+ Main props are id, interval and n_intervals
		* id
		* interval -> time in ms between refreshes
		* n_intervals -> counts the number of refreshes done, increasing +1 for each. Very very useful for the callback function (as you pass the prop n_intervals as input, and every time it refreshes, the callback is called)
			- ex. `@app.callback(Output('storage', 'data'), Input('timer', 'n_intervals'))`
		* max_intervals -> maximum number of refreshes allowed




## Dash Bootstrap Components (DBC)
`import dash_bootstrap_components as dbc` -> requires pip installation

DBC is a package that helps with the layout of the app, giving pre-styled bootstrap components. It's useful but in the long run it causes problems, I'd rather not use it and create a custom CSS for the company to use

This package has several components similar to the dash HTML components, but they work better between them. Pretty much pick one or the other, you're going to have compatibility issues if you mix.







--------------------------------------------------



# Dash Callbacks
Callbacks allow user interactivity within the app. It's the mechanism that connects the Dash components to each other (eg once user selects dropdown value, thje figure is updated).
Callbacks are composed of two parts: decorator (identifies the relevant components for the action) and callback function (defines how dash components should interact)
The callback structure can be analysed with the development tools on the webapp page, the blue '<>' button. This is only available in debug mode (debug=True), remember to remove when deploying!!




## Callback decorator
The callback decorator registers the callback function, telling it when to call the function and where to send the return value of the function. 
It basically is an event listener that redirects the input to the callback function and takes the output of the callback function to redirect it to the desired component. Ex. dropdown element gets selected, callback function makes graph based on it, decorator redirects graph data from callback function to correct dcc.Graph

The callback decorator is placed above callback function, with no space between the two. The decorator takes two arguments, _Output and Input_. Input is 'where' the decorator should listen for a change, and output is where the output of the callback function should be redirected.
_Input and Output_ thake in two arguments:
- _component_id_ -> corresponds to the ID of a particular component
- _component_property_ -> corresponds to the specific prop of that component that has to be targetted 

There can be more than 1 inputs and outputs
ex. Select dropdown element to produce a new chart 
```py
@app.callback(
	Output(component_id='chart', component_property='figure'), # figure is the prop that stores the graph itself
	[Input(component_id='dropdown', component_property='value')]
	)
```


#### State in callback decorator - Multiple inputs
When creating a complex graph, this might require several inputs from different components. For example the following callback:
```py
@app.callback( 
	Output("my-choropleth", "figure"), 
	Input("my-button", "n_clicks"), 
	Input("storage", "data"), 
	State("years-range", "value"), 
	State("radio-indicator", "value"),
)
```
The callback outputs a graph in the 'figure' prop of 'my-choropleth'; to do this it requires the input of the button, or the refresh of the data in the storage. In order to capture other data, such as the range of years to display, or which values from the radio selection to include, we use the 'State' instead of 'Input'. The 'State' doesn't trigger the callback when its components are altered, but keeps track of the user's selection. When the callback is then triggered by one of the inputs, the callback function can have access to the value of the two other components represented in the state!




## Callback function
The callback function is what handles the 'event' registered by the callback decorator. It takes in as many inputs as given, and returns as many objects as required, following the order of input/return.
(eg if there are 2 outputs in the callback, the function's return statement will return value 1 to the first output and value 2 to the second. ex. dropdown value changes, the callback picks up on it and sends the new value of the dropdown to the callback fn, the callback produces the graph and sends it to the component in the first output, it then produces a string describing the graph and sends it to the textbox underneath the graph in the second output)
(The same is true for inputs, if there are 3 inputs in a callback, you need to interact with all 3 to trigger the callback)
Usually these functions handle null values and return null strings/lists/dicts accordingly.
No need for an example as it's literally a function. _Just be careful about_ not altering global variables, you can duplicate them by filtering DFs or deepcopying them with the module deepcopy.



## Block callback triggering
By default, all callbacks are triggered when the app starts. To stop this behaviour (for example for callbacks that produce graphs, even if the DF is not loaded), there are 2 ways:
- prevent_initial_callbacks=True -> in the Dash() declaration (applies to all I suppose)
	+ app = Dash(`__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], prevent_initial_callbacks=True`)
- prevent_initial_call=True -> in specific callback (applies only to the callback where it's added)
	+ @app.callback(`Output("storage", "data"), Input("timer", "n_intervals"), prevent_initial_call=True`)













---------------------------

# Plotly Express
Plotly express is a high level interface for creation of graphs. This is the 'easier' version of plotly, with, as trade-off, a lower level of customisability.

Plotly requires the data to be formatted with pandas' DataFrames, or in a format interpretable by pandas. For ease every time there's written DF it refers to a Pandas' DataFrame (or data convertible to it), and DF columns refer to the strings composing the DF's columns names OR array/Series-like objects which can hold the same function


### Line graphs
Line graphs are among the easiest graphs. The three main props required are 'data_frame', 'x' and 'y' (specify which columns plot and in which axis). 
`plotly.express.line(data_frame=df, x='axis_name', y='axis_name'`

These are all the props for this graph, some of them might be useful, I'll do a better breakdown of what each one controls later on!:

```py
plotly.express.line(data_frame=None, x=None, y=None, 
	title=None, # Title of the figure
	template=None, #Figure template name, search more in plotly site
	labels={}, # Rename DF columns to be displayed with better form. Takes in a dict in format {'col_name':'New label'}
	text=None, # Column name. Values from this column appear in the figure as text labels
	error_x=None, error_y=None, # Column name. Used to size x/y-axis error bars in the positive direction
	error_x_minus=None, error_y_minus=None,  # Column name. Used to size x/y-axis error bars in the negative direction
	orientation=None, # String. 'h' for horizontal, 'v' for vertical (def 'v')
	color=None, # Categorical column name. Differentiate groups of data with colour based on category
	color_discrete_sequence=None, # List of strings. Strings fo valid CSS colours. When prop color is set, colour values are assigned by cycling through the ones listed in this prop (if category orders is set it'll follow that order)
	color_discrete_map={}, # Dictionary {'values_of_color(prop)_column': 'color value'}. Assign directly colours to the categories
	line_dash=None, # Column name. Values are used to assign dash-patterns to lines
	line_dash_sequence=None, # List of strings. Strings defining valud plotly dash patterns. When line_dash is set, the dash-patterns are assigned by cycling though this list (if category_orders is set, it'll follow that order) 
	line_dash_map={}, # Dictionary {'values_of_line_dash(prop)_column': dash_pattern}. Map dash patterns directly to the categorical variable.
	range_x=None, range_y=None, # List of two numbers. Overrides auto scaling and sets min and max value to be displayed
	line_shape=None, # String. Either 'linear'(default) or 'spline' (?) 
	render_mode='auto', # String. Either 'auto' (default), 'svg', 'webgl'. 'svg' for imgs with <1000 data points. 'webgl' >1000 data points. 'auto' picks between one of the other two modes
	category_orders={}, # Used to forcefully order categorical values in columns. Takes in dictionary in form of {'col_name': ['list', 'of', 'values', 'in', 'desired', 'order']}
	width=None, height=None, # max dimensional values for graph in px
	log_x=False, log_y=False, # Present axis in a logarithmic scale (log10 by default can be personalised with more complex graph objects)
	hover_name=None, # Column name. Values in column will appear in bold in the hover tooltip
	hover_data=None, # List of column names or dictionary with column names as keys and True or formatting string (ex 3f) or list like data. Dictates what extra data will be presented on hover. Example hover_data=['column_name']
	custom_data=None, # Names of columns or array-like objects. These are extra meta data not visible to the user, but included in events emitted by the figure 
	line_group=None, # Column name. Column values are used to group rows of data into lines (??)
	facet_row=None, facet_col=None, # Column name. Values are used to assign marks to faccetted subplots in the vertical/horizontal direction respectively
	facet_col_wrap=0, # Integer. Maximum number of facet columns. Ignored if 0
	facet_row_spacing=None, facet_col_spacing=None, # Float 0<x<1. Spacing between facet rows, in paper units (???). Defaults, row:0.03, col:0.02
	animation_frame=None, # Column name. Values are used to assign marks to animation frames
	animation_group=None # Column name. Values are used to provide constancy across animation frames; rows with same animation_group are treated as if they describe the same object (???)
	)



```



























































































Callback diagram page 125














































