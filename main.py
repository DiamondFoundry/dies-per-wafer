import functions_framework
import google.cloud.logging
import logging

logging_client = google.cloud.logging.Client()
logging_client.setup_logging()


### CLOUD FUNCTIONS ENTRYPOINT
@functions_framework.http
def cloudFunctionEntrypoint(request):
    '''
    Entry point for Google Cloud Functions
        This function must be chosen as the entrypoint when defining the function in the Google Cloud Console
    '''
    s = 'Hello world'
    logging.debug(s)
    return s
