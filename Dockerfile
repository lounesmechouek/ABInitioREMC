FROM python:3.10.11

# set the working directory
WORKDIR /app

# copy the config file to the working directory
COPY ./server_config/config.toml /root/.streamlit/config.toml

# run server configuration
#RUN python config.py

# copy the files to the container
COPY . /app

# install dependencies
RUN pip install .

# Expose port 8501 for Streamlit
EXPOSE 8501

# Start the app
CMD ["streamlit", "run", "run.py"]