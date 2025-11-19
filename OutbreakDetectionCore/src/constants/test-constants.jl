export CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC,
    CLINICAL_TEST_SPECS, PROTOTYPE_RDT_TEST_SPECS

const CLINICAL_CASE_TEST_SPEC = IndividualTestSpecification(1.0, 0.0, 0)
const EPI_LINKED_CASE_TEST_SPEC = IndividualTestSpecification(1.0, 0.8, 0)
const CLINICAL_TEST_SPECS = (CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC)
const PROTOTYPE_RDT_TEST_SPECS = (
    IndividualTestSpecification(0.43, 0.98, 0),
    IndividualTestSpecification(0.26, 0.97, 0),
    IndividualTestSpecification(0.95, 0.95, 0),
    IndividualTestSpecification(0.75, 0.96, 0),
)
