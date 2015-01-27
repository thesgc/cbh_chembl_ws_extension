Feature: CBH Compound Batch
    In order to keep track of compounds I would like to create 
    a batch of a new compound in the project I am an admin of.
    As a project stakeholder I do not want batches to be logged in 
    projects that the user does not have editor access to
#behave && cat typescript | aha > test.html

    Scenario Outline: Batches project privileges saving and validating
        Given a User
        and I have a valid molfile
        and I have a valid list of SMILES
        and I have an invalid Excel File
        and a valid project exists proja
        and I automatically have editor permissions as creator
        and I remove my permissions
        and I have<editor> given myself editor privileges for proja
        and I have<viewer> given myself viewer privileges for proja        
        When <dologin> log in
        Then I <action> my cbh_compound_batch to proja and the <responsecode>


        Examples: Validation api
        |   dologin    |   editor  |   viewer  |   responsecode   | action |
        |      I         |           |           |  response code will be 202             | validate |
        |       I        |           |   nt     |  response code will be 202             | validate |
        |        I       |   nt     |           |  response code will be 401             | validate |
        |         I      |   nt     |   nt     |  response code will be 401             | validate |
        |   I do not      |           |           |  response code will be 401             | validate |
        
        Examples: Create api
             |   dologin    |   editor  |   viewer  |   responsecode   | action |
         |      I         |           |           |  response code will be 201            | create |
        |       I        |           |   nt     |  response code will be 201             | create |
        |        I       |   nt     |           |  response code will be 401             | create |
        |         I      |   nt     |   nt     |  response code will be 401             | create |
        |   I do not      |           |           |  response code will be 401             | create |




    Scenario Outline: Batches project privileges get list
        Given a User
        and I have a valid molfile
        and I have a valid list of SMILES
        and I have an invalid Excel File
        and a valid project exists proja
        and I automatically have editor permissions as creator
        and a single batch exists in proja
        and I remove my permissions
        and I have<editor> given myself editor privileges for proja
        and I have<viewer> given myself viewer privileges for proja        
        When <dologin> log in
        Then I <action> my cbh_compound_batch to proja and the <responsecode>


        
        Examples: Get List api  - Note the viewer can list
        |   dologin    |   editor  |   viewer  |   responsecode   | action |
         |      I         |           |           |   response code will be 1_memberlist      | list |
        |       I        |           |   nt     |   response code will be 1_memberlist      | list |
        |        I       |   nt     |           |   response code will be 1_memberlist   | list |   
        |         I      |   nt     |   nt     |   response code will be 0_memberlist      | list |

        
        Examples: Get api  - Note the viewer can list
        |   dologin    |   editor  |   viewer  |   responsecode   | action |
         |      I         |           |           |   response code will be 200     | get |
        |       I        |           |   nt     |   response code will be 200      | get |
        |        I       |   nt     |           |   response code will be 200   | get |   
        |         I      |   nt     |   nt     |   response code will be 401      | get |















        #Examples: Validate multi batch



        #Examples: Upload File - to be done manually





    Scenario: User submits a substance to project they have editor rights to, but the substance is already preregistered to other private and public projects they do not have editor rights to.  
        Given a User
        and I have a valid molfile
        and a valid project exists proja
        and a valid project exists projb
        and a valid project exists projc
        and a valid project exists projd
        and this substance is not pre-registered to proja 
        and this substance is pre-registered to projb as private
        and this substance is pre-registered to projc as private
        and this substance is pre-registered to projd as public
        and I have editor rights for proja
        and I have viewer rights for projb
        and I have no rights for projc
        and I have no rights for projd
        When I log in 
        and I validate a cbh_compound_batch to proja        
        Then the response will contain an id for a previously registered substance in projb
        and the response will contain an id for a previously registered substance in projd
        and the response will not contain an id for a previously registered substance in projc


    Scenario: User submits retrieved ID of substance to make a batch, but is not entitled to this substance.
        Given a User
        When I log in 
        and I have a valid molfile
        and a valid project exists proja
        and a valid project exists projb
        and a valid project exists projc
        and a valid project exists projd
        and this substance is not pre-registered to proja 
        and this substance is pre-registered to projb as private
        and this substance is pre-registered to projc as private
        and this substance is pre-registered to projd as public
        and I have editor rights for proja
        and I have viewer rights for projb
        and I have no rights for projc
        and I have no rights for projd
        When I validate a cbh_compound_batch to proja        
        Then the response will contain an id for a previously registered substance in projb
        and the response will not contain an id for a previously registered substance in projc
        and the response will contain an id for a previously registered substance in projd
        When I register substance as a batch to retrieved ID
        Then the response for projb registered ID will be <responsecode>
        and the response for projd registered ID will be <responsecode>
       

    



    Scenario: User registers a batch with a substance ID that does not exist 
        Given a User
        When I log in 
        and I have a substance ID <false_ID> that does not exist
        When I create batch with <false_ID>       
        Then the response will be <error_substance_ID_not_found>




    Scenario: User registers a batch with a substance ID that is not valid 
        Given a User
        When I log in 
        and I have a substance ID <invalid_ID> that is not valid
        When I create batch with <invalid_ID>       
        Then the response will be <error_invalid_ID_recieved>












    """
    Created using
    arrays = [("","do not"), ("","not"), ("","not"), ("401",)]
    rows = list(itertools.product(*arrays))
    rows = ["|\t" +  "\t\t|\t".join(row) + "\t\t|" for row in   list(itertools.product(*arrays))]
    print "\n".join(rows)
    """





    #Scenario Outline: Other project privileges, ivladifd project key

