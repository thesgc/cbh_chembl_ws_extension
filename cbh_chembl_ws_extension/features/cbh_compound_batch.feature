Feature: CBH Compound Batch
    In order to keep track of compounds I would like to create 
    a batch of a new compound in the project I am an admin of.
    As a project stakeholder I do not want batches to be logged in 
    projects that the user does not have editor access to


    Scenario Outline: Batches project privileges
        Given a User
        When <dologin> log in
        and I have a valid molfile
        and a valid project exists proja
        and I have <editor> editor privileges for proja
        and I have <viewer> viewer privileges for proja
        when I save my cbh_compound_batch to proja
        Then the response code will be <responsecode>

        Examples: Tests
        |   dologin    |   editor  |   viewer  |   responsecode   | 
        |      I         |           |           |   201             |
        |       I        |           |   not     |   201             |
        |        I       |   not     |           |   401             |
        |         I      |   not     |   not     |   401             |
        |   do not      |           |           |   401             |
 


    Scenario: User submits a substance to project they have editor rights to, but the substance is already preregistered to other private and public projects they do not have editor rights to.  
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
        when I validate a cbh_compound_batch to proja        
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



    Scenario:  User registers a L-chiral_substance when D-chiral_subtance and rmix_chiral_substance already exist.
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














    """
    Created using
    arrays = [("","do not"), ("","not"), ("","not"), ("401",)]
    rows = list(itertools.product(*arrays))
    rows = ["|\t" +  "\t\t|\t".join(row) + "\t\t|" for row in   list(itertools.product(*arrays))]
    print "\n".join(rows)
    """





    #Scenario Outline: Other project privileges, ivladifd project key

