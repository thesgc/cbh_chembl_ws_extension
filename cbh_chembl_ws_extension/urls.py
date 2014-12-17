from django.conf.urls import patterns, url
from cbh_chembl_ws_extension import api as webservices
from cbh_chembl_ws_extension import Login
from chembl_webservices import __version__ as ws_version
from chembl_webservices import api_name
from chembl_core_db.utils import DirectTemplateView
from django.conf import settings


spore_context={
    "WS_VERSION": ws_version,
    "WS_BASE_URL": settings.WS_BASE_URL,
    "WS_DOCS_TITLE": settings.WS_DOCS_TITLE
}

urlpatterns = patterns('',
    url(r'chemblws/login^$',Login.as_view(), name="login"),
    url(r'^%s/docs' % api_name, DirectTemplateView.as_view(template_name="docs.html"), name='ws_docs'),
    url(r'^%s/spore' % api_name, DirectTemplateView.as_view(template_name="ws_spore.json" , extra_context=spore_context), name='ws_spore_endpoint'),
)
urlpatterns += webservices.urls