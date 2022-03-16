from django.shortcuts import render
from django.http import HttpResponse
from .forms import UploadFileForm
import sys
from django.http import HttpResponseRedirect
#from .ion_predictor.django_wrapper.auto_predictor import screen_predict
import os

# Create your views here.


def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")


"""

# ------------------------------------------------------------------
def file_upload(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            handle_uploaded_file(request.FILES['file'])
            file_obj = request.FILES['file']
            path="media/documents/"+file_obj.name

            #calculate 
            df=screen_predict(path)

            #remove file
            try:
                os.remove(path)
            except:
                pass

            #return csv
            response = HttpResponse(content_type='text/csv; charset=utf8')
            response['Content-Disposition'] = 'attachment; filename=users.csv'
            df.to_csv(path_or_buf=response, encoding='utf_8_sig', index=None)
            return response 

    else:
        form = UploadFileForm()
    return render(request, 'ion_manager/upload.html', {'form': form})
#
#
# ------------------------------------------------------------------
def handle_uploaded_file(file_obj):
    file_path = 'media/documents/' + file_obj.name 
    with open(file_path, 'wb+') as destination:
        for chunk in file_obj.chunks():
            destination.write(chunk)
"""
