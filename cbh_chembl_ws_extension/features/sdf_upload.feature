SDF upload	

Feature: SDF upload
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
        The the responce will be <This substance has already been registered. Would you like to register as a new batch (or force registration)?>


	








        