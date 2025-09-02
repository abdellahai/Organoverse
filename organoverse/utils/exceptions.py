"""
Custom exceptions for Organoverse workbench.
"""


class OrganoverseError(Exception):
    """Base exception for Organoverse workbench."""
    pass


class ModuleError(OrganoverseError):
    """Exception raised by individual modules."""
    pass


class ValidationError(OrganoverseError):
    """Exception raised during input validation."""
    pass


class ConfigurationError(OrganoverseError):
    """Exception raised for configuration issues."""
    pass


class DatabaseError(OrganoverseError):
    """Exception raised for database access issues."""
    pass


class AssemblyError(OrganoverseError):
    """Exception raised during assembly process."""
    pass


class ModelError(OrganoverseError):
    """Exception raised by AI/ML models."""
    pass