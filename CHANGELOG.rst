^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package bayes_tracking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1.0.2 (2014-09-03)
------------------
* Adding boost as a build an run dependency.
* Contributors: Christian Dondrup

1.0.1 (2014-09-02)
------------------
* Merge pull request `#4 <https://github.com/cdondrup/bayestracking/issues/4>`_ from cdondrup/master
  Renamed `bayestracking` to `bayes_tracking` to comply with ros naming conventions.
* Renamed `bayestracking` to `bayes_tracking` to comply with ros naming conventions.
* Merge pull request `#2 <https://github.com/cdondrup/bayestracking/issues/2>`_ from cdondrup/master
  Adding optional catkin build
* * Enabling normal cmake build if catkin is not available.
  * Updated package.xml and create install targets for ROS release
  * Refactored README to markdown
* Refactored bayestracking library to be a catkin package
  Created a cleaner file structure by seperating the src from the include files.
* Update README
* Update README
* Now generates cmake config files to enable usage of find_package
* Added the creation of a pkg-config file to make the library easier to use.
* PFilter not working right now in simple_tracking.cpp, but all other tested filters (UKF and EKF) worked very well
* fixed size_t uint issue
* compiles now (but simple_tracking example still gives errors)
* initial import form package
* Contributors: Christian Dondrup, Marc Hanheide
