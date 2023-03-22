import time

def mock_output(files, buffer):
    for file in files:
        print('Starting file ', file)#, file=buffer)
        buffer.append(f'Starting file {file}')
        time.sleep(1)
        print('Complete')#, file=buffer)
        



