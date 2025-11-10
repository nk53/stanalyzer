import argparse
import collections.abc as abc
import sys
import typing as t
from functools import wraps
from pathlib import Path

import pydantic
from pydantic import BaseModel, ValidationError
from pydantic.fields import FieldInfo
from pydantic_core._pydantic_core import PydanticUndefined, PydanticUndefinedType

from stanalyzer.cli import cli_input
from stanalyzer.validation import Project
from stanalyzer.utils import write_settings

# Generics
BM = t.TypeVar('BM', bound=BaseModel)
IM = t.TypeVar('IM', bound='InteractiveModel')

# Pydantic field definition types: sequences
FieldList: t.TypeAlias = abc.Sequence[str] | None
OptFieldList: t.TypeAlias = FieldList | None

# Pydantic field definition types: mappings
FieldDict: t.TypeAlias = abc.MutableMapping[str, t.Any]
FieldInfoDict: t.TypeAlias = abc.MutableMapping[str, t.Any]

# InteractiveModel internal types
FactoryNoArgs: t.TypeAlias = abc.Callable[[], t.Any]
FactoryArgs: t.TypeAlias = abc.Callable[[dict[str, t.Any]], t.Any]
Factory: t.TypeAlias = FactoryNoArgs | FactoryArgs
OptFactory: t.TypeAlias = FactoryNoArgs | None

# InteractiveModel arguments
ModelInject: t.TypeAlias = type[BaseModel] | abc.Sequence[type[BaseModel]] | None


def undefined(value: object) -> t.TypeGuard[PydanticUndefinedType | None]:
    if value is PydanticUndefined or value is None:
        return True
    return False


class InteractiveModel(BaseModel):
    """BaseModel that interactively validates missing input values from stdin.

    Example usage:
        from datetime import date
        from pydantic import Field

        class MyModel(InteractiveModel):
            id: int = 0
            name: str
            crated_date: date = Field(default_factory=date.today)


        no_interaction_needed = MyModel(id=1, name='Normal BaseModel',
                        created_date=date(2000,1,1))  # same as BaseModel
        some_interaction = MyModel(name='Some Interaction')
        # example interaction session (after ': ' is user supplied input):
        #   id (0): 2
        #   created_date (auto):
        print(*list(some_interaction))  # expected output below:
        # ('id', 2) ('name', 'Some Interaction') ('created_date', date(2000, 1, 1))
    """
    # must remain unset, else __init_subclass__ overrides subclass attrs
    _sta_field_models: t.ClassVar[list[type[BaseModel]]]
    _sta_no_write_fields: t.ClassVar[list[str]]
    _sta_no_interact_fields: t.ClassVar[list[str]]
    _sta_default_overrides: t.ClassVar[dict[str, object]]

    @classmethod
    def as_model(cls: type[IM], other: type[BM], *args: t.Any, **kwargs: t.Any) -> BM:
        """Shortcut to initialize interactively and return non-interactive obj"""
        while True:
            try:
                obj = cls(*args, **kwargs)
                model_dict = obj.model_dump()
                return other(**model_dict)
            except ValidationError as exc:
                errors = exc.errors()
                print("Got errors:")
                for error in errors:
                    locs = ', '.join(map(str, error['loc']))
                    print(f"    {locs}: {error['msg']}")

                    # remove failed inputs
                    for loc in error['loc']:
                        if not isinstance(loc, str):
                            raise NotImplementedError(repr(loc))
                        model_dict.pop(loc, None)

                # retry asking only for failed inputs
                kwargs = model_dict
            except (KeyboardInterrupt, EOFError):
                sys.exit(1)

    def __init__(self, *args: t.Any, **kwargs: t.Any):
        super().__init__(*args, **kwargs)

    @classmethod
    def _sta_setdefault(cls) -> None:
        for attr, ann in InteractiveModel.__annotations__.items():
            value = getattr(cls, attr, PydanticUndefined)

            if isinstance(value, PydanticUndefinedType):
                value = getattr(cls.__private_attributes__, attr, PydanticUndefined)
            if isinstance(value, PydanticUndefinedType):
                ann_cls = ann.__args__[0]
                setattr(cls, attr, ann_cls())
            if isinstance(value, pydantic.fields.ModelPrivateAttr):
                setattr(cls, attr, value.default)

    def __init_subclass__(cls, _sta_inject: ModelInject = None, **kwargs: t.Any):
        """Use _sta_inject to inherit options from a BaseModel.

        This can be used to override the validation behavior of an existing BaseModel
        without the need to redefine it as a subclass of InteractiveModel.

        Parameters
        ==========
            cls          InteractiveModel descendant
            _sta_inject  a BaseModel or list of BaseModels to inherit options from;
                         later classes' options override earlier ones

        Parameter specification
        =======================
        """
        cls._sta_setdefault()
        if _sta_inject is None:
            _sta_inject = tuple()
        elif isinstance(_sta_inject, abc.Sequence):
            _sta_inject = tuple(_sta_inject)
        elif issubclass(_sta_inject, BaseModel):
            _sta_inject = (_sta_inject,)

        value: PydanticUndefinedType | FieldInfo | t.Any
        for _inject in _sta_inject:
            for field, ann in getattr(_inject, '__annotations__', {}).items():
                if field in cls._sta_no_write_fields:
                    continue
                value = _inject.model_fields.get(field, PydanticUndefined)
                cls.__annotations__[field] = ann
                if isinstance(value, FieldInfo):
                    model_field = value
                else:
                    model_field = FieldInfo.from_annotated_attribute(ann, value)

                if field in cls._sta_default_overrides:
                    default = cls._sta_default_overrides[field]
                    if callable(default):
                        model_field = FieldInfo.merge_field_infos(
                            model_field, default=PydanticUndefined,
                            default_factory=default)
                    else:
                        model_field = FieldInfo.merge_field_infos(
                            model_field, default=default,
                            default_factory=None)

                cls.model_fields[field] = model_field
                setattr(cls, field, model_field)

        def default_factory(params: FieldDict, other_factory: OptFactory = None) -> Factory:
            def get_result() -> object:
                result = cli_input.get_input(**params)
                if type(result) in cls._sta_field_models:
                    fields = list(result)
                    if len(fields) > 1:
                        raise NotImplementedError("field dependencies")

                    return fields[0][1]

                if result == '' and other_factory is not None:
                    return other_factory()

                return result

            return get_result

        for field, ann in getattr(cls, '__annotations__', {}).items():
            if field in cls._sta_no_interact_fields:
                continue  # let pydantic handle it

            value = getattr(cls, field, PydanticUndefined)
            params: FieldDict = {'prompt': f"{field}: "}
            default_str: object = PydanticUndefined
            default_func: OptFactory = None

            if isinstance(value, PydanticUndefinedType):
                value = cls.model_fields.get(field, PydanticUndefined)

            if isinstance(value, PydanticUndefinedType):
                info = FieldInfo.from_annotation(ann)
            elif isinstance(value, FieldInfo):
                # just add the annotation
                value.annotation = ann
                info = value

                if callable(info.default_factory):
                    default_str = 'auto'
                    params['default'] = ''
                    default_func = t.cast(FactoryNoArgs, info.default_factory)
                elif info.default is not PydanticUndefined:
                    default_str = info.default  # allows default=None
                    params['default'] = info.default
            else:
                default_str = value
                params['default'] = value
                info = FieldInfo.from_annotated_attribute(ann, value)

            if default_str is not PydanticUndefined:
                if default_str == '':
                    default_str = repr('')
                params['prompt'] = f"{field} ({default_str}): "

            params['dtype'] = cls._sta_wrap__({field: info})

            new_field = pydantic.Field(
                default_factory=default_factory(params, default_func))

            new_field.annotation = info.annotation
            setattr(cls, field, new_field)

        super().__init_subclass__(**kwargs)

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs: t.Any) -> None:
        # TODO: cls.__pydantic_decorators__.field_validators
        #       cls.__pydantic_decorators__.model_validators
        super().__pydantic_init_subclass__(**kwargs)

    @classmethod
    def _sta_wrap__(cls: type[IM], fields: FieldDict) -> abc.Callable[..., IM]:
        field_defs: FieldInfoDict = {
            name: (info.annotation, info)
            for name, info in fields.items()
        }
        # field_model = pydantic.create_model(  # type: ignore[call-overload]
        field_model = pydantic.create_model(
            f'{"_".join(fields.keys())}_model', **field_defs)

        @wraps(cls)
        def wrapper(*args, **pydantic_kwargs):
            """Accepts pydantic way of initializing model or positional-only"""
            if pydantic_kwargs:
                return field_model(**pydantic_kwargs)

            return field_model(**{k: v for k, v in zip(fields, args)})

        cls._sta_field_models.append(field_model)
        return wrapper


class InteractiveProject(InteractiveModel, _sta_inject=Project):
    _sta_no_interact_fields = ['id', 'SLURM', 'PBS']
    _sta_no_write_fields = ['id']
    _sta_default_overrides = {
        'output_path': Path('.').resolve(),
        'input_path': Path('.').resolve(),
    }

    @classmethod
    def write_settings(cls: type['InteractiveProject'], output_path: str = 'project.json',
                       other: type[Project] = Project, *args: t.Any, **kwargs: t.Any) -> None:
        # NB: see TODO in InteractiveModel.__pydantic_init_subclass__
        model = cls.as_model(other)
        exclude = set(cls._sta_no_write_fields)
        write_settings(
            Path(model.input_path) / output_path,
            model.model_dump(mode='json', exclude=exclude))


def main(output_path: str = 'project.json') -> None:
    parser = argparse.ArgumentParser(
        prog="stanalyzer config",
        description="Interactively create project.json")
    parser.parse_args(sys.argv[2:])

    InteractiveProject.write_settings(output_path)
