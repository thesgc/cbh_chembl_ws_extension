from behave import given, when, then
import json
from django.contrib.auth.models import User, Group


@given('I have a valid molfile')
def step(context):
    context.post_data["ctab"] = """


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


@given("I have a valid list of SMILES")
def step(context):
    pass


@given("I have an invalid Excel File")
def step(context):
    pass


@given("a single batch exists in {projkey}")
def step(context, projkey):
    from cbh_core_model.models import Project
    from cbh_chembl_model_extension.models import CBHCompoundBatch
    p = Project.objects.get(project_key=projkey)
    context.batch = CBHCompoundBatch.objects.create(
        ctab=context.post_data["ctab"], project=p)


@then("I {action} my cbh_compound_batch to {projkey} and the response code will be {responsecode}")
def step(context, action=None, projkey=None, responsecode=None):
    from cbh_core_model.models import Project
    if action in ["validate", "create"]:
        if action == "validate":
            path = "/dev/cbh_compound_batches/validate/"
            func = context.api_client.post
        elif action == "create":
            path = "/dev/cbh_compound_batches/"
            func = context.api_client.post
        context.post_data["projectKey"] = projkey
        resp = func(
            path,
            format='json',
            data=context.post_data,
        )
        print(resp.status_code)
        print(resp.__dict__)
        assert resp.status_code == int(responsecode)
    else:
        from cbh_core_model.models import Project

        path = "/dev/cbh_compound_batches/"
        if action == "get":
            path = path + str(context.batch.id)
            # print(path)
        resp = context.api_client.get(
            path,
        )
        print(resp.status_code)
        print(resp.__dict__)
        if action == "get":
            assert resp.status_code == int(responsecode)
        else:
            project_id = Project.objects.get(project_key=projkey).id

            assert len(context.ser.deserialize(resp.content)["objects"]) == int(
                responsecode.split("_")[0])
