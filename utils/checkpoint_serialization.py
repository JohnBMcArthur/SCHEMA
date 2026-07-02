"""
Serialization utilities for checkpoint system.
Handles complex data structures that can't be directly pickled or JSON-serialized.
"""

import json
import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Union

try:
    import schemarecomb
    from schemarecomb import Library
    SCHEMARECOMB_AVAILABLE = True
except ImportError:
    SCHEMARECOMB_AVAILABLE = False
    Library = None


def serialize_library(library) -> str:
    """
    Serialize a schemarecomb Library object to JSON string.
    
    Args:
        library: schemarecomb.Library object
        
    Returns:
        str: JSON string representation
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required for Library serialization")
    
    if hasattr(library, 'to_json'):
        return library.to_json()
    else:
        # Fallback: convert to dict
        return json.dumps({
            'breakpoints': [bp.position for bp in library.breakpoints],
            'energy': str(library.energy),
            'mutation_rate': str(library.mutation_rate),
            'gg_prob': str(library.gg_prob) if hasattr(library, 'gg_prob') else None,
        })


def deserialize_library(json_str: str) -> Any:
    """
    Deserialize a schemarecomb Library object from JSON string.
    
    Args:
        json_str: JSON string representation
        
    Returns:
        schemarecomb.Library object
    """
    if not SCHEMARECOMB_AVAILABLE:
        raise ImportError("schemarecomb is required for Library deserialization")
    
    try:
        return Library.from_json(json_str)
    except Exception as e:
        # Try fallback deserialization
        data = json.loads(json_str)
        raise ValueError(f"Could not deserialize Library: {str(e)}")


def serialize_library_list(libraries: List[Any]) -> List[str]:
    """
    Serialize a list of Library objects.
    
    Args:
        libraries: List of Library objects
        
    Returns:
        List of JSON strings
    """
    if not libraries:
        return []
    
    result = []
    for lib in libraries:
        if lib is None:
            result.append(None)
        elif SCHEMARECOMB_AVAILABLE and isinstance(lib, Library):
            result.append(serialize_library(lib))
        else:
            # Try to serialize as-is (might be tuple or other format)
            result.append(json.dumps(lib, default=str))
    return result


def deserialize_library_list(json_list: List[str]) -> List[Any]:
    """
    Deserialize a list of Library objects from JSON strings.
    
    Args:
        json_list: List of JSON strings
        
    Returns:
        List of Library objects
    """
    if not json_list:
        return []
    
    result = []
    for json_str in json_list:
        if json_str is None:
            result.append(None)
        else:
            try:
                result.append(deserialize_library(json_str))
            except Exception as e:
                import logging
                logging.getLogger(__name__).warning("deserialize_library_list: deserialize_library failed: %s", e)
                try:
                    result.append(json.loads(json_str))
                except Exception as e2:
                    logging.getLogger(__name__).warning("deserialize_library_list: json fallback failed: %s", e2)
                    result.append(None)
    return result


def serialize_numpy_array(arr: np.ndarray) -> Dict[str, Any]:
    """
    Convert numpy array to JSON-serializable dict.
    
    Args:
        arr: numpy array
        
    Returns:
        dict: Serializable representation
    """
    return {
        '__numpy_array__': True,
        'data': arr.tolist(),
        'dtype': str(arr.dtype),
        'shape': arr.shape
    }


def deserialize_numpy_array(data: Dict[str, Any]) -> np.ndarray:
    """
    Reconstruct numpy array from serialized dict.
    
    Args:
        data: Serialized array dict
        
    Returns:
        numpy array
    """
    return np.array(data['data'], dtype=data['dtype']).reshape(data['shape'])


def make_json_safe(obj: Any) -> Any:
    """
    Convert object to JSON-safe format.
    Handles numpy arrays, sets, and other non-JSON types.
    
    Args:
        obj: Object to convert
        
    Returns:
        JSON-safe representation
    """
    if isinstance(obj, np.ndarray):
        return serialize_numpy_array(obj)
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif isinstance(obj, np.generic):
        return obj.item()
    elif isinstance(obj, set):
        return {'__set__': True, 'items': list(obj)}
    elif isinstance(obj, (tuple, list)):
        return [make_json_safe(item) for item in obj]
    elif isinstance(obj, dict):
        return {key: make_json_safe(value) for key, value in obj.items()}
    elif hasattr(obj, '__dict__'):
        # Try to serialize object as dict
        return {'__object__': True, 'class': obj.__class__.__name__, 'data': make_json_safe(obj.__dict__)}
    else:
        return obj


def restore_from_json_safe(obj: Any) -> Any:
    """
    Restore object from JSON-safe format.
    
    Args:
        obj: JSON-safe representation
        
    Returns:
        Restored object
    """
    if isinstance(obj, dict):
        if '__numpy_array__' in obj:
            return deserialize_numpy_array(obj)
        elif '__set__' in obj:
            return set(obj['items'])
        elif '__object__' in obj:
            # Can't fully restore arbitrary objects, return dict
            return obj['data']
        else:
            return {key: restore_from_json_safe(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [restore_from_json_safe(item) for item in obj]
    else:
        return obj


def serialize_raspp_raw_results(raw_results: Any) -> Any:
    """
    Serialize RASPP raw results (may be Library objects or tuples).
    
    Args:
        raw_results: Raw RASPP results (list of Library objects or tuples)
        
    Returns:
        Serialized representation
    """
    if raw_results is None:
        return None
    
    if isinstance(raw_results, list):
        # Check if it's a list of Library objects
        if raw_results and SCHEMARECOMB_AVAILABLE and isinstance(raw_results[0], Library):
            return {
                '__library_list__': True,
                'libraries': serialize_library_list(raw_results)
            }
        else:
            # Try to serialize as tuples or other formats
            return [make_json_safe(item) for item in raw_results]
    else:
        return make_json_safe(raw_results)


def deserialize_raspp_raw_results(serialized: Any) -> Any:
    """
    Deserialize RASPP raw results.
    
    Args:
        serialized: Serialized representation
        
    Returns:
        Original format (Library objects or tuples)
    """
    if serialized is None:
        return None
    
    if isinstance(serialized, dict) and serialized.get('__library_list__'):
        return deserialize_library_list(serialized['libraries'])
    elif isinstance(serialized, list):
        return [restore_from_json_safe(item) for item in serialized]
    else:
        return restore_from_json_safe(serialized)
