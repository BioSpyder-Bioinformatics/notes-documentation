import tkinter as tk
import tkinter.ttk as tkk
import tkinter.font as tkFont
from tkinter.filedialog import askopenfilenames, askdirectory
from mock_functions import mock_output
import threading
import io
import os
import time

# For guide:
# - All files have to have the same format! (eg fastq or fastq.gz)


########################
# Functions
def get_reference_list():
    #this_dir = os.getcwd()
    # Get path of this file
    this_file = os.path.realpath(__file__)
    # From file path get relative directory
    this_dir = this_file.removesuffix(os.path.basename(this_file)) # removes this filename from path
    files = os.listdir(f'{this_dir}references/')
    files = [file.removesuffix('.fa').removeprefix('TempO-Seq_') for file in files if file.split('.')[-1] == 'fa']
    return files

def get_path_from_reference(reference):
    # Get path of this file
    this_file = os.path.realpath(__file__)
    # From file path get relative directory
    this_dir = this_file.removesuffix('aligner_gui.py')
    return_path = f'{this_dir}references/TempO-Seq_{reference}.fa'
    return return_path


def update_output(process, buffer, files):
    # Number of completed files
    total_completed = 0
    total_files = len(files)
    # Update output while the process is alive (in separate fn bc need to thread)
    while process.is_alive():
        print('alive', process.is_alive())  
        # If the buffer is longer than one, consume it and append it to output of interest
        while len(buffer) > 0:
            # Pop first element and append it to output of interest
            text = buffer.pop(0)
            # Update progress bar if keyword (this case is starting) in buffer
            if 'Starting' in text:
                total_completed += 1
                progress_bar7['value'] = (total_completed/total_files)*100

            # Enable editing of text6
            text6['state'] = 'normal'
            # Set new text6 (1.0 is starting row 1, column 0; up to tk.END (end of text); text to set)
            text6.insert(tk.END, text+'\n')
            # Disable editing
            text6['state'] = 'disabled'
        print(buffer)  
        time.sleep(1)


# Main function controlling everything
def submit_btn():
    # Retrieve variables 
    # Reference genome
    reference_genome = variable1.get()
    # Get actual path to genome
    reference_genome = get_path_from_reference(reference_genome)
    
    # Aligner
    aligner = variable2.get()
    # Make aligner lowercase!
    aligner = aligner.lower()

    # Number of mismatches
    mismatches = variable3.get()

    # Number of threads
    threads = variable4.get()
    
    # Get files from files output (in string format)
    files_text = text5.get('1.0', tk.END)
    # Format strings
    files = files_text.split('\n')
    files = [file.strip() for file in files if file != ''] # strip spaces and remove empty files
    print(files)
    
    # Set communicating output box (This probably is best removing and replacing )
    # output_label = tkk.Label(frame3, text='Empty')
    # # Assign it
    # output_label.grid(row=4)

    # If there are no files, let the user know that you need at least one!!
    if len(files) < 1 :
        text6['state'] = 'normal'
        text6.replace('1.0', tk.END, 'Please select at least one file to start!')
        text6['state'] = 'disabled'
        return
    else:
        # Clean up text6 
        text6['state'] = 'normal'
        text6.replace('1.0', tk.END, '')
        text6['state'] = 'disabled'



    # Declare string buffer to capture stdout
    buffer = [] #io.StringIO()

    # This is so the window does not freeze
    #mock_output(files)
    # Before sending the function, determine if files are zipped or not!
    process = threading.Thread(target=mock_output, args=[files, buffer]) # Need to give a buffer to append output to
    # This in a try catch?? catch displays error message!
    process.start()

    # Also so the window does not freeze, updates the output separately
    threading.Thread(target=update_output, args=[process, buffer, files]).start()

    print('hello', reference_genome, aligner)



def get_directory_or_files():
    # Connected to button, prompts the user to get the files
    var = askopenfilenames() or askdirectory() 
    # Format text to display nicely
    text = ' \n'.join([x for x in var])
    # Enable editing of text5
    text5['state'] = 'normal'
    # Set new text5 (1.0 is starting row 1, column 0; up to tk.END (end of text); text to set)
    text5.replace('1.0', tk.END, text)
    # Disable editing
    text5['state'] = 'disabled'
    






########################


# 3 Quadrants (frames)
# q1 -> All dropdowns + checkbox stacked up
# q2 -> select files + progress output
# q3 -> (horizontal) progress bar + start btn

# Standard initialisation
# Create a new window
window = tk.Tk()
# Set the window title
window.title("TempoSeq Aligner")
# Set the window size
window.geometry("1200x800") #width x height


# Set up font objects
arial20 = tkFont.Font(family='Arial', size=20)
arial15 = tkFont.Font(family='Arial', size=15)


# Declare frames and assign them to grid!
frame1 = tk.LabelFrame(window, height=680, width=600, bg='red', padx=100, pady=200) # Play around with padding n stuff
frame1.grid(row=0, column=0)
# frame1.configure(padding=5, borderwidth=2, relief='solid')

frame2 = tk.LabelFrame(window, height=680, width=600, bg='blue')
frame2.grid(row=0, column=1)

frame3 = tk.LabelFrame(window, height=150, width=1200, bg='purple')
frame3.grid(row=1, column=0, columnspan=2)


# Quadrant 1 (column 0) -----------------------------------
# Elements by #
# 1 -> select reference genome #
# 2 -> select aligner
# 3 -> select n mismatches
# 4 -> select n threads 

# Reference genome label + dropdown
#genome_options = ['rat_w_1.0', 'human_s1500_1.2', 'human_w_2.0', 'human_w_2.1', 'mouse_s1500_1.2', 'mouse_w_1.0']
genome_options = get_reference_list()
label1 = tkk.Label(frame1, text='Select reference genome', font=arial20)
label1.grid(row=0)
# reference genome variable
variable1 = tk.StringVar(frame1, genome_options[0]) 
dropdown1 = tk.OptionMenu(frame1, variable1, *genome_options)
dropdown1.config(font=arial15, justify='center', width=35)
dropdown1.grid(row=1)

# Select aligner label + dropdown
aligner_options = ['STAR', 'BWA', 'Kallisto']
label2 = tkk.Label(frame1, text='Select aligner', font=arial20)
label2.grid(row=2)
# Selected aligner variable
variable2 = tk.StringVar(frame1, aligner_options[0])
dropdown2 = tk.OptionMenu(frame1, variable2, *aligner_options)
dropdown2.config(font=arial15, justify='center', width=20)
dropdown2.grid(row=3)

# Select number mismatches label + dropdown (star only)
mismatches_options = ['0', '1', '2', '3']
label3 = tkk.Label(frame1, text='Select number of mismatches (STAR only)', font=arial20)
label3.grid(row=4)
# Mismatches variable
variable3 = tk.StringVar(frame1, mismatches_options[2])
dropdown3 = tk.OptionMenu(frame1, variable3, *mismatches_options)
dropdown3.config(font=arial15, justify='center', width=20)
dropdown3.grid(row=5)

# Select number of threads label + spinbox
label4 = tkk.Label(frame1, text='Select the number of threads', font=arial20)
label4.grid(row=6)
# Threads variable
variable4 = tk.StringVar(frame1, value=8)
spinbox4 = tkk.Spinbox(frame1, from_=1, to=50, textvariable=variable4)
spinbox4.grid(row=7)


#----------------------------------------------------------

# Quadrant 2 (column 1) -----------------------------------
# Elements by # 
# 5 -> select files 
# 6 -> progress output#

# Select files label + button + textbox 
label5 = tkk.Label(frame2, text='Select directory or files', font=arial20)
label5.grid(row=0, column=0)
button5 = tkk.Button(frame2, text='Select', command=get_directory_or_files)
button5.grid(row=0, column=1)
# Textbox variable to be updated
text5 = tk.Text(frame2)#, width=50, height=50)
# Don't allow input!
text5['state'] = 'disabled'

text5.grid(row=1, column=0, columnspan=2)


# Progress output label + textbox
label6 = tkk.Label(frame2, text='Processing output: ', font=arial20)
label6.grid(row=2, column=0, columnspan=2)
text6 = tk.Text(frame2)
text6['state'] = 'disabled'
text6.grid(row=3, column=0, columnspan=2)

#----------------------------------------------------------

# Quadrant 3 (row 1) -----------------------------------
# Elements by #
# 7 -> progress bar
# 8 -> run button

# Progress bar (spans like 4/5ths of the bottom)
progress_bar7 = tkk.Progressbar(frame3, length=500, maximum=100, mode='determinate')
progress_bar7.grid(column=0, columnspan=5)

# Start button 
button8 = tkk.Button(frame3, text='Start alignment', command=submit_btn)
button8.grid(column=6)



if __name__ == '__main__':
    # Start loop 
    window.mainloop()





# Need to do:
# File upload
# Progress bar!
