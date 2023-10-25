# Use an official Ubuntu as a base image
FROM ubuntu:latest

# Set environment variables to non-interactive (this avoids some prompts)
ENV DEBIAN_FRONTEND noninteractive

# Install Python and pip
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    python3-tk\
    unixodbc-dev \
    g++ 

# Install Microsoft ODBC Driver for SQL Server
RUN apt-get install -y wget apt-transport-https
RUN wget https://packages.microsoft.com/keys/microsoft.asc -O- | apt-key add -
RUN wget https://packages.microsoft.com/config/ubuntu/$(. /etc/os-release; echo $VERSION_ID)/prod.list
RUN mv prod.list /etc/apt/sources.list.d/mssql-release.list
RUN apt-get update
RUN ACCEPT_EULA=Y apt-get install -y msodbcsql17

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Install any needed packages specified in requirements.txt
RUN pip3 install --no-cache-dir -r requirements.txt

# Make port 80 available to the world outside this container
EXPOSE 5025

# Define environment variable
ENV NAME World



CMD ["python", "runserver.py"]