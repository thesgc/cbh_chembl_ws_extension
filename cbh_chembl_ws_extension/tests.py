#Tests to implement

'''
ModelResources will be created for each of the required models

Therefore we will need to creaste tests for:

CBHCompoundResource - will need to check for permissions via batches

MoleculeBatchResource - which will also update compounds and molecules and also create and update permissions on batchs

Project??

GroupResource

UserResource

'''

import datetime
from django.contrib.auth.models import User
from tastypie.test import ResourceTestCase
from cbh_chembl_model_extension.models import CBHCompoundBatch
from django.db import connection


class CompoundBatchResourceTest(ResourceTestCase):
    # Use ``fixtures`` & ``urls`` as normal. See Django's ``TestCase``
    # documentation for the gory details.
    #fixtures = ['test_entries.json']

    def setUp(self):
        super(CompoundBatchResourceTest, self).setUp()
        # Create a user.
        self.username = 'daniel'
        self.password = 'pass'
        self.user = User.objects.create_user(self.username, 'daniel@example.com', self.password)

        # Fetch the ``Entry`` object we'll use in testing.
        # Note that we aren't using PKs because they can change depending
        # on what other tests are running.

        # We also build a detail URI, since we will be using it all over.
        # DRY, baby. DRY.
        # self.detail_url = '/api/v1/entry/{0}/'.format(self.entry_1.pk)

        # The data we'll send on POST requests. Again, because we'll use it
        # frequently (enough).
        self.post_data = {
            'ctab' : """


  8  8  0  0  0  0            999 V2000
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0
    0.8660    0.5000    0.0000 C   0  0  0  0  0  0
    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0
    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0
   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0
    0.0000    2.0000    0.0000 O   0  0  0  0  0  0
    0.0000   -2.0000    0.0000 O   0  0  0  0  0  0
  1  21  0     0  0
  2  32  0     0  0
  3  41  0     0  0
  4  51  0     0  0
  5  62  0     0  0
  6  11  0     0  0
  1  72  0     0  0
  4  82  0     0  0
M  END"""

        ,
            'editable_by': {"fdsfsdf":4324234},
            'viewable_by': {},
        }

    def setup_session(self):
        self.api_client.client.login(username=self.username, password=self.password)

    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)

    def test_post_list_unauthenticated(self):
        self.assertHttpUnauthorized(self.api_client.post('/chemblws/cbh_compound_batches/', format='json', data=self.post_data))

    def test_post_list(self):
        self.setup_session()

        # Check how many are there first.
        self.assertEqual(CBHCompoundBatch.objects.count(), 0)
        self.assertHttpCreated(self.api_client.post('/chemblws/cbh_compound_batches/', format='json', data=self.post_data))
        # Verify a new one has been added.
        self.assertEqual(CBHCompoundBatch.objects.count(), 1)



    def test_post_list_validate(self):
        self.setup_session()

        # Check how many are there first.
        resp = self.api_client.post('/chemblws/cbh_compound_batches/validate/', format='json', data=self.post_data)
        self.assertHttpAccepted(resp)


    def test_post_list(self):
        self.setup_session()

        # Check how many are there first.
        self.assertEqual(CBHCompoundBatch.objects.count(),0)
        resp = self.api_client.post('/chemblws/cbh_compound_batches/', format='json', data=self.post_data)
        self.assertHttpCreated(resp)
        self.assertEqual(CBHCompoundBatch.objects.count(),1)
        c = CBHCompoundBatch.objects.filter(warnings__contains={"pains_count":"1"})
        self.assertEqual(c.count(),1)
        print resp.__dict__

