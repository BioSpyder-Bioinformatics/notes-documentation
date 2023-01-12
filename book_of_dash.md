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
- Variables and DFs declared outside a function are global!
- The layout of the app refers to the alignment of the components. The style refers to how the element look, also known as *props*
- Dash can render _cytoscape_, a js library for rendering complex graphs. It can also render VTK for 3D graphs


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



## Callback function
The callback function is what handles the 'event' registered by the callback decorator. It takes in as many inputs as given, and returns as many objects as required, following the order of input/return.
(eg if there are 2 outputs in the callback, the function's return statement will return value 1 to the first output and value 2 to the second. ex. dropdown value changes, the callback picks up on it and sends the new value of the dropdown to the callback fn, the callback produces the graph and sends it to the component in the first output, it then produces a string describing the graph and sends it to the textbox underneath the graph in the second output)
(The same is true for inputs, if there are 3 inputs in a callback, you need to interact with all 3 to trigger the callback)
Usually these functions handle null values and return null strings/lists/dicts accordingly.
No need for an example as it's literally a function. _Just be careful about_ not altering global variables, you can duplicate them by filtering DFs or deepcopying them with the module deepcopy.





---------------------------

# Plotly Express















































































































































