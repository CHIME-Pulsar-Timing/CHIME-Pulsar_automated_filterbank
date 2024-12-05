import sys
import logging

# LOGGING UNHANDLED EXCEPTIONS:
# from https://stackoverflow.com/a/16993115/18740127
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.critical(
        "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
    )
