from django.conf.urls import patterns, url, include
from cbh_chembl_ws_extension import api as webservices
from cbh_chembl_ws_extension import Login, Logout
from chembl_webservices import __version__ as ws_version
from chembl_webservices import api_name
from chembl_core_db.utils import DirectTemplateView
from django.conf import settings
from flowjs import urls as flow
from django.contrib import admin

admin.autodiscover()


urlpatterns = patterns('',
    url(r'^%s/' % api_name ,Login.as_view(), name="login1"),
    url(r'^%s/login' % api_name ,Login.as_view(), name="login"),
    url(r'^%s/logout' % api_name ,Logout.as_view(), name="logout"),
    url(r'^%s/flow/' % api_name, include(flow)), #adding this to allow configured upload URL within django-flowjs
    url(r'^%s/admin/' % api_name, include(admin.site.urls)),
	url(r'^grappelli/', include('grappelli.urls')),
)
urlpatterns += webservices.urls


if "django_webauth" in settings.INSTALLED_APPS:
    urlpatterns += patterns('',
                            url(r'^%s/webauth/' % api_name, include('django_webauth.urls', 'webauth')),
                            )