FROM python:3

WORKDIR /usr/src/app

#do requirements before code because they change less
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .
RUN pip install --no-cache-dir .

# TODO change this
ENTRYPOINT ["TODO"]
# TODO change this
CMD ["TODO"]
