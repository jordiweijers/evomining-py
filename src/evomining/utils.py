import logging

def setup_logging(log_file, name="main_logger", console=True):
    """ 
    Setup up logging.
    Args:
        log_file (str): The path to the log file.
        name (str): The name of the logger.
        console (bool): Whether to also log to the console.
    Returns:
        logging.Logger: The configured logger.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    if logger.hasHandlers():
        logger.handlers.clear()
    if console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter("%(message)s")
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
	)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    return logger