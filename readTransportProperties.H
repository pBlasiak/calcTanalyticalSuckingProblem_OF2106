IOdictionary SuckingProblemDict
(
	IOobject
	(
		"SuckingProblemDict",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// Phase change type
word phaseChangeType(SuckingProblemDict.lookup("phaseChangeType"));
// Automatic set of starting value for epsilon
bool epsilonStartingValue(SuckingProblemDict.lookupOrDefault("epsilonStartingValue", false));
// Superheating temperature
dimensionedScalar Tinf("Tinf", dimTemperature, SuckingProblemDict);
// Wall temperature
dimensionedScalar Tw("Tw", dimTemperature, SuckingProblemDict);
