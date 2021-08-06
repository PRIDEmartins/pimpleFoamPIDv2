#include "createRDeltaT.H"

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Info << "Reading field pidEp\n" << endl;
volScalarField pidEp
(
    IOobject
    (
        "pidEp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define user-defined PID constants that will be read by the solver     //
// from the file ./constant/pidProperties                                //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading pidProperties\n" << endl;

IOdictionary pidProperties
(
    IOobject
    (
        "pidProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar pidKp // Proportional error gain
(
    "pidKp",
    dimensionSet( 0, 0, 0, 0, 0, 0, 0 ),
    pidProperties.lookup("pidKp")
);
Info<< "Proportional error gain (Force coefficient)" << endl;


dimensionedScalar pidTi // PID Integration time
(
    "pidTi",
    dimensionSet( 0, 0, 1, 0, 0, 0, 0 ),
    pidProperties.lookup("pidTi")
);
Info<< "PID Integration time" << endl;


dimensionedScalar pidTd // PID Derivative time
(
    "pidTd",
    dimensionSet( 0, 0, 1, 0, 0, 0, 0 ),
    pidProperties.lookup("pidTd")
);
Info<< "PID Derivative time" << endl;


dimensionedScalar pidNum // Number of assimilation points: nudging points in the system
(
    "pidNum",
    dimensionSet( 0, 0, 0, 0, 0, 0, 0 ),
    pidProperties.lookup("pidNum")
);
Info<< "Number of assimilation points\n" << endl;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define user-defined PID constants that will be read by the solver     //
// from the file ./constant/pidProperties                                //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Info<< "Reading pidData\n" << endl;

IOdictionary pidData
(
    IOobject
    (
        "pidData",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<vector> pidPositions
(
    pidData.lookup("pidPositions") // location of the nudging points in the mesh
);

List<scalar> pidTime
(
    pidData.lookup("pidTime") // time at which the data was collected
);

List<scalar> pidP0
(
    pidData.lookup("pidP0")  // reference pressure (e.g., at inlet) at the time-step t
);

List<scalar> pidP
(
    pidData.lookup("pidP") // Pressure values at the pidNum - 1 locations at which pressure will be nudged. Values shold be places such that: piP(point1(0:nt);point2(0:nt);...)
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define other quantities
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Locate closest point in the grid relative to pidPositions
std::vector<int> pidPositions_inMesh(pidNum.value());
for (int i = 0; i < pidNum.value(); ++i) {
    pidPositions_inMesh[i] = mesh.findCell(pidPositions[i]);
}

// Define the pid time vector as the input ./constants/pidData->pidTime
std::vector<double> pidTime_vector(pidTime.size());
for (int i = 0; i < pidTime.size(); ++i) {
    pidTime_vector[i] = pidTime[i];
}

// Convert the input ./constants/pidData->pidP0 into a matrix that constains as i-th column the reference pressure value, pidP0, at the t-th time-step: pidP_matrix(length(pidTime),length(pidNum-1)
// Firstly, create a matrix of size  length(pidNum-1) x length(pidTime)
std::vector< std::vector<double> > pidP_matrix(pidNum.value()-1);
for ( int i = 0 ; i < (pidNum.value()-1) ; i++ ) {
    pidP_matrix[i].resize( pidP0.size() );
}

// Secondly, the i-th row represents the i-th nudging point, and each column represents the t-th time step
for (int i = 0; i < (pidNum.value()-1); ++i) {
    for (int j = 0; j < ( pidP0.size() ); ++j) {
        pidP_matrix[i][j] = pidP[ j + i*( pidP0.size() ) ];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Display info about points where the pressure is being corrected
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "PID: point pressure data\n" << endl;
Info<< "Reference point position = " << pidPositions[0] << nl << endl;

for (int i=1; i < (pidNum.value()); ++i) {
    Info<< "Assimilation point " << i << endl;
    Info<< "Position = " << pidPositions[i] << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

// Copied from OFv6.0:
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"