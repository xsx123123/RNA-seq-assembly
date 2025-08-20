import sys
from loguru import logger
from datetime import datetime

logger.remove()
log_file_name = f'{config['project_name']}_{datetime.now().strftime("%Y-%m-%d-%H:%M:%S")}.log'
new_level = logger.level("SNAKY", no=38, color="<yellow>", icon="🐍")
logger.add(log_file_name,rotation="500 MB",
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> <blue> LOG-LEVELS: </blue>"
                  "{level.icon} <blue> LOG-INFO : </blue> "
                  "{message}",
           level=config["log_level"],
           colorize=True)
# logger example usage
# logger.info("This is an info message.")
# logger.debug("This is a debug message.")
# logger.warning("This is a warning message.")
# logger.error("This is an error message.")
# logger.success("This is a success message.")