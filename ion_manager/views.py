from django.shortcuts import render
from django.http import HttpResponse
from .forms import UploadFileForm
import os

# Create your views here.


def index(request):
    return HttpResponse("Hello, world.")


