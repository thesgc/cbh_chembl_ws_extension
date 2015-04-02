
Feature: SDF upload
'''
    I wish to submit an SDF file to a project I have editor rights. The substances in the SDF file are either pre-registered or not registered publically on the CBH resistration system.

	Scenario: User submits an SDF file containing unregistered substances to proja. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and a substance in the SDF file is not pre-registered
        and I have editor rights for proja
        when I submit file
        then the SDF substances will be registered to proja

    Scenario: User submits an SDF file containing pre-registered substances to proja. 
        Given a User
        When I log in 
        and a valid project exists proja
        and a substance in the SDF file is pre-registered
        and I have editor rights for proja
        when I submit file
        Then the response will be This substance has already been registered. Would you like to register as a new batch (or force registration)

    Scenario: User submits a SDF file containing unregistered substances to proja. Some of the substances contain compounds with one stereocentre.
        Given a User
        When I log in
        and a valid project exists proja
        and I have editor rights for proja
        and a substance in the SDF file contains one stereocentre
        when I submit file
        then the substances will be registered to proja and retain their single stereocentre

    Scenario: User submits a SDF file containing unregistered substances to proja. Some of the substances contain compounds with multiple stereocentres.
        Given a User
        When I log in
        and a valid project exists proja
        and I have editor rights for proja
        and a substance in the SDF file contains multiple stereocentres
        when I submit file
        then the substances will be registered to proja and retain their multiple stereocentres

    Scenario: User submits a SDF file containing pre-registered substances to proja. Some of the pre-registered substances contain compounds with one stereocentre.
        Given a User
        When I log in
        and a valid project exists proja
        and a substance containing one stereocentre in the SDF file is pre-registered with matching stereochemistry
        and I have editor rights for proja
        and a substance in the SDF file contains stereochemistry
        when I submit file
        Then the response will be This substance has already been registered. Would you like to register as a new batch (or force registration)

    Scenario: User submits a SDF file containing pre-registered substances to proja. Some of the pre-registered substances contain compounds with multiple stereocentres.
        Given a User
        When I log in
        and a valid project exists proja
        and a substance containing multiple stereocentres in the SDF file is pre-registered with matching stereochemistry
        and I have editor rights for proja
        and a substance in the SDF file contains stereochemistry
        when I submit file
        Then the response will be This substance has already been registered. Would you like to register as a new batch (or force registration)







	



'''




        