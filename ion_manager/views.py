from django.shortcuts import render
from django.http import HttpResponse
from .forms import UploadFileForm
from django.http import HttpResponse
import sys
from django.http import HttpResponseRedirect
from .ion_predictor.django_wrapper.auto_predictor import screen_predict

# Create your views here.
def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")





# ------------------------------------------------------------------
def file_upload(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            sys.stderr.write("*** file_upload *** aaa ***\n")
            handle_uploaded_file(request.FILES['file'])
            file_obj = request.FILES['file']
            sys.stderr.write(file_obj.name + "\n")
            df=screen_predict(file_obj.name)

            response = HttpResponse(content_type='text/csv; charset=utf8')
            response['Content-Disposition'] = 'attachment; filename=users.csv'
            df.to_csv(path_or_buf=response, encoding='utf_8_sig', index=None)
            return response 

            return HttpResponseRedirect('predict')
    else:
        form = UploadFileForm()
    return render(request, 'ion_manager/upload.html', {'form': form})
#
#
# ------------------------------------------------------------------
def handle_uploaded_file(file_obj):
    sys.stderr.write("*** handle_uploaded_file *** aaa ***\n")
    sys.stderr.write(file_obj.name + "\n")
    file_path = 'media/documents/' + file_obj.name 
    sys.stderr.write(file_path + "\n")
    with open(file_path, 'wb+') as destination:
        for chunk in file_obj.chunks():
            sys.stderr.write("*** handle_uploaded_file *** ccc ***\n")
            destination.write(chunk)
            sys.stderr.write("*** handle_uploaded_file *** eee ***\n")
#
# ------------------------------------------------------------------
def predict(request):
    str_out = "Success!<p />"
    str_out += "成功<p />"
    return HttpResponse(str_out)
# ------------------------------------------------------------------