import sys
import requests

topic = "mudock-lorenzo-8f3a92c7e1"
message = sys.argv[1] if len(sys.argv) > 1 else "Experiment finished"

requests.post(
    f"https://ntfy.sh/{topic}",
    data=message.encode("utf-8"),
)