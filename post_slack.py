from slack_sdk import WebClient
import os
import glob
pngs = glob.glob("nsub*/*.png")

#store token as an environment variable
token = os.environ['SLACK_API_TOKEN']
client = WebClient(token=token)
for im in pngs:
    response = client.files_upload(channels="#chimepulsar-lpt", file=im, title="Candidate", initial_comment="Candidate")
