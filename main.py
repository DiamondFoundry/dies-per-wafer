import functions_framework
import google.cloud.logging
import logging
from opt.dies_per_wafer_latest import DiesPerWaferCalculator

logging_client = google.cloud.logging.Client()
logging_client.setup_logging()


def fit_solution(request_json):
    try:
        width = request_json["width"]
        height = request_json["height"]
        xspacing = request_json["xspacing"]
        yspacing = request_json["yspacing"]
        waferdiameter = request_json["input_wafer_diameter"]
        edgeexclusionwidth = request_json["edge_exclusion_width"]
        ft_grid = request_json["ft_grid"]
        symmetric = request_json["symmetric"]
        ft_ShiftRows = request_json["ft_ShiftRows"]
        ft_ShiftCols = request_json["ft_ShiftCols"]
        ft_ShiftRot = request_json["ft_ShiftRot"]
    except Exception as e:
        logging.error(f"Missing entry in request: {e}")
        raise e

    dpw = DiesPerWaferCalculator(
        width=width,
        height=height,
        xspacing=xspacing,
        yspacing=yspacing,
        waferdiameter=waferdiameter,
        edgeexclusionwidth=edgeexclusionwidth,
        ft_grid=ft_grid,
        searchdepth=0,
        symmetric=symmetric,
        ft_ShiftRows=ft_ShiftRows,
        ft_ShiftCols=ft_ShiftCols,
        ft_ShiftRot=ft_ShiftRot,
    )
    dpw.fit()
    json_str = dpw.format_json_obj()
    return json_str


### CLOUD FUNCTIONS ENTRYPOINT
@functions_framework.http
def calculate_dies_per_wafer(request):
    """
    Entry point for Google Cloud Functions
        This function must be chosen as the entrypoint when defining the function in the Google Cloud Console
    """
    request_json = request.get_json()
    result = fit_solution(request_json)
    return result


if __name__ == "__main__":
    test_input = {
        "width": 20.0,
        "height": 20.0,
        "xspacing": 1.0,
        "yspacing": 1.0,
        "edge_exclusion_width": 1.0,
        "input_wafer_diameter": 200.0,
        "ft_ShiftRows": True,
        "ft_ShiftCols": True,
        "ft_ShiftRot": False,
        "symmetric": False,
        "ft_grid": True,
    }
    result = fit_solution(test_input)
    print(result)
