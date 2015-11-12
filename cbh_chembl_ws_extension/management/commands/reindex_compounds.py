from django.core.management.base import BaseCommand, CommandError
from django.http import HttpRequest

class Command(BaseCommand):

    def handle(self, *args, **options):
        from cbh_chembl_ws_extension.compounds import CBHCompoundBatchResource

        cbr = CBHCompoundBatchResource()
        cbr.reindex_elasticsearch(HttpRequest())
