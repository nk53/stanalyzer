import sys
import typing as t
from pydantic import ValidationError

__all__ = ['get_choice', 'get_input', 'validate_type']


def get_bool(prompt: str, default: bool = False, retry_prompt: t.Optional[str] = None):
    if default:
        prompt = f"{prompt} [Y/n] "
        default_str = 'y'
    else:
        prompt = f"{prompt} [y/N] "
        default_str = 'n'

    allowed = ['y', 'ye', 'yes', 'n', 'no']
    yes_no_response = get_input(prompt, default=default_str, case_sensitive=False,
                                retry_prompt=retry_prompt, allowed=allowed)

    if yes_no_response.lower().startswith('y'):
        return True

    return False


def get_choice(prompt: str, choices: t.Sequence, descriptions: t.Optional[t.Sequence[str]] = None,
               dtype: type = str, default: t.Any = None, retry_prompt: t.Optional[str] = None):
    choices = [str(choice) for choice in choices]
    tpl = "[{}]{}"
    if descriptions is None:
        descriptions = [''] * len(choices)
    else:
        descriptions = list(' '+s for s in descriptions)

    for choice, description in zip(choices, descriptions):
        print(tpl.format(choice, description))
    return get_input(prompt, dtype=dtype, allowed=choices,
                     retry_prompt=retry_prompt, default=default)


def get_input(prompt: str, dtype: type = str, default: t.Any = None,
              allowed: t.Optional[t.Sequence] = None, ntries: t.Optional[int] = None,
              case_sensitive: bool = True, retry_prompt: t.Optional[str] = None):
    if retry_prompt is None:
        retry_prompt = "Invalid response, please try again, or interrupt to quit: "

    if not case_sensitive and allowed:
        allowed = [s.lower() for s in allowed]

    validate_type(retry_prompt, str, 'retry_prompt')
    if ntries is not None:
        validate_type(ntries, int, 'ntries')

    while ntries is None or ntries > 0:
        try:
            response = input(prompt)
            if response == '':
                if default is not None:
                    return default

            return dtype(response)

            if allowed is None:
                return response
            elif response in allowed:
                return response
            elif not case_sensitive and response.lower() in allowed:
                return response
        except (EOFError, KeyboardInterrupt):
            print()
            sys.exit(1)
        except ValidationError as exc:
            # get first error message
            errors = exc.errors()
            error = errors[0]
            if len(error['loc']) == 1:
                if ctx := error.get('ctx'):
                    errmsg = ctx.get('error', error.get('msg'))
                else:
                    errmsg = error['msg']
                print(errmsg)
                prompt = "Please try again, or interrupt to quit: "
            elif len(error['loc']) > 1:
                for error in errors:
                    if len(loc := error['loc']) > 2:
                        raise NotImplementedError(repr(loc))

                    cls, attr = loc[:2]
                    msg = error['msg']
                    print(f"{attr}: {msg}")

                if isinstance(cls, int):
                    raise NotImplementedError(repr(loc))

                if info := getattr(dtype, cls, None):
                    if examples := info.examples:
                        examples = ', '.join(repr(ex) for ex in examples)
                        print(f"{cls} examples: {examples}")
            else:
                raise Exception("No loc data")
        except Exception as exc:
            prompt = retry_prompt
            print("Got unexpected exception:", repr(exc))
        finally:
            if ntries is not None:
                ntries -= 1
    print("Too many failed attempts. Exiting")
    sys.exit(2)


def validate_type(obj: t.Any, cls: type, var_name: str):
    tpl = "Invalid type for {}; expected {}, got '{}'"
    if not isinstance(obj, cls):
        expected = type(cls).__name__
        got = type(obj).__name__
        raise TypeError(tpl.format(var_name, expected, got))
