/*
 * File:   testEvergreen.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on Aug 4, 2018, 4:09:20 PM
 */

#ifndef TESTANEN_H
#define TESTANEN_H

#include <cppunit/extensions/HelperMacros.h>

class testEvergreen : public CPPUNIT_NS::TestFixture {

    CPPUNIT_TEST_SUITE(testEvergreen);
    CPPUNIT_TEST(testSoltrack_);
    CPPUNIT_TEST(testPvwattsv5_);
    CPPUNIT_TEST_SUITE_END();

public:
    testEvergreen();
    virtual ~testEvergreen();

private:
    void testSoltrack_();
    void testPvwattsv5_();
};

#endif /* TESTANEN_H */
