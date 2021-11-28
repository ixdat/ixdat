# Test structure

The tests are organized into folders, where each folder is for a specific type of test
(such as unit tests or functional tests). The filenames (when possible) refer to the module
that contains the unit under test. In this way, there may be a identically named file in
several folders containing different types of tests for the same module, as shown below.

```
tests/
|- functional/
   |- __init__.py
   |- test_measurements.py
|- unit/
   |- __init__.py
   |- test_measurements.py
```
