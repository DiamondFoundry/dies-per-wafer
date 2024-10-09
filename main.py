import functions_framework
import google.cloud.logging
import logging
import json

logging_client = google.cloud.logging.Client()
logging_client.setup_logging()

def main(request_json): 
    id_user = request_json.get('id_user', None)
    s = 'Hello world'
    logging.debug(s)
    return json.dumps({'content': {'s': s, 'id_user': id_user}})


### CLOUD FUNCTIONS ENTRYPOINT
@functions_framework.http
def cloudFunctionEntrypoint(request):
    '''
    Entry point for Google Cloud Functions
        This function must be chosen as the entrypoint when defining the function in the Google Cloud Console
    '''
    request_json = request.get_json()
    result = main(request_json)
    return result

if __name__ == '__main__':
    test_input = {"id_user": 549}
    result = main(test_input)
    print(result)
