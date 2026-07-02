"""
Simple test script to verify imports work correctly.
"""

import sys
from pathlib import Path

print("Testing imports...")

try:
    print("1. Testing schema_raspp package...")
    from schema_raspp import schema, raspp, pdb
    print("   ✓ schema_raspp imports successful")
except Exception as e:
    print(f"   ✗ Error importing schema_raspp: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

try:
    print("2. Testing utils modules...")
    from utils import file_handlers, schema_wrapper, raspp_wrapper, visualization
    print("   ✓ utils imports successful")
except Exception as e:
    print(f"   ✗ Error importing utils: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

try:
    print("3. Testing Streamlit...")
    import streamlit as st
    print("   ✓ streamlit import successful")
except Exception as e:
    print(f"   ✗ Error importing streamlit: {e}")
    print("   Note: This is OK if streamlit is not installed in this environment")

try:
    print("4. Testing other dependencies...")
    import numpy
    import matplotlib
    import plotly
    import pandas
    print("   ✓ All dependencies imported successfully")
except Exception as e:
    print(f"   ✗ Error importing dependencies: {e}")
    import traceback
    traceback.print_exc()

print("\n✓ Import tests completed!")
