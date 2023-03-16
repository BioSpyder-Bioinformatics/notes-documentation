#import win32com.client as win32 # conda install -c conda-forge pywin32 #pip install -i https://pypi.anaconda.org/rida/simple pypiwin32

def send_email(recipient):
    import smtplib

    body = 'Subject: Alignment complete \n\n' + 'Hello,\n\nYour alignment has completed!' + '\n\nRegards,\nBioSpyder Team'
    sender = 'biospyder_aligner@outlook.com'

    # Initialise smtp objects
    try:
        smtpObj = smtplib.SMTP('smtp-mail.outlook.com', 587)  #smtp.office365.com
    except Exception as e:
        print(e)
        smtpObj = smtplib.SMTP_SSL('smtp-mail.outlook.com', 465)

    smtpObj.ehlo()
    smtpObj.starttls()
    smtpObj.login('biospyder_aligner@outlook.com', 'Biospyder123_')
    smtpObj.sendmail(sender, recipient, body)
    smtpObj.quit()
    print(f'Email sent to {recipient}')



recipient =  'salvocamiolo@biospyder.com' #gioelemook97@gmail.com' #'gioelemontis@biospyder.com' #'salvocamiolo@biospyer.com' 

send_email(recipient)

