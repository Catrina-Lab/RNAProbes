# A collection of utility methods and classes specifically for this project

import zipfile
from argparse import Namespace
from io import UnsupportedOperation
from pathlib import Path
import io

from src import util
from src.util import remove_if_exists

class UnclosableStringIO(io.StringIO):
    def close(self):
        # Override close so it doesn't actually close the stream
        pass

class UnclosableBytesIO(io.BytesIO):
    def close(self):
        # Override close so it doesn't actually close the stream
        pass


class ProgramObject:
    def __init__(self, output_dir: Path, file_stem: str, arguments: Namespace, **kwargs):
        self.result_obj = Namespace(**kwargs)
        self.output_dir = output_dir
        self.file_stem = file_stem
        self.arguments = arguments
        if output_dir is not None: output_dir.mkdir(parents=True, exist_ok=True)

    def save_buffer(self, rel_path: str):
        """
        Returns a buffer to use to store files
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :return:
        """
        return self.file_path(rel_path) #test

    def open_buffer(self, rel_path: str, mode = "w"):
        return open(self.save_buffer(rel_path), mode)

    def reset_buffer(self, rel_path: str):
        remove_if_exists(self.file_path(rel_path))

    def file_path(self, rel_path: str) -> Path:
        """
        Returns a path to use to store files. Same as save_buffer if from the command line, but will always
        return a path (unlike save_buffer which may return a buffer)
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :return:
        """
        return self.output_dir / self._format_relative_path(rel_path)

    def _format_relative_path(self, rel_path: str) -> str:
        return rel_path.replace("[fname]", self.file_stem)
    def register_file(self, rel_path: str, true_path:  Path = None):
        pass
    def get_arg(self, argument):
        return getattr(self.arguments, argument)

    def set_args(self, **kwargs):
        for argument, value in kwargs.items():
            setattr(self.arguments, argument, value)

    def get_result_arg(self, argument):
        return getattr(self.result_obj, argument)

    def set_result_args(self, **kwargs):
        for argument, value in kwargs.items():
            setattr(self.result_obj, argument, value)
    def to_zip(self, name) -> tuple[bytes, str]:
        archive_name = name.replace("[fname]", self.file_stem)
        return util.get_folder_as_zip(self.output_dir), archive_name

#todo: prevent collision in file_dict and buffer_dict
class BufferedProgramObject(ProgramObject):
    def __init__(self, output_dir: Path, file_stem: str, arguments: Namespace, **kwargs):
        ProgramObject.__init__(self, output_dir, file_stem, arguments, **kwargs)
        self.buffer_dict = dict()
        self.file_dict = dict() #separate so as not to mess with other methods using buffer_dict
    def save_buffer(self, rel_path: str, is_string: bool = True):
        """
        Returns a string buffer to use to store files
        :param rel_path: the relative path to store files in. Replaces [fname] with the file stem
        :param is_string: True if using StringIO, false for BytesIO
        :return:
        """
        path = Path(self._format_relative_path(rel_path))
        return self.buffer_dict.setdefault(path, UnclosableStringIO() if is_string else UnclosableBytesIO()) #test

    def open_buffer(self, rel_path: str, mode = "w"):
        """

        :param rel_path:
        :param mode: the mode to write in, w to remove the previous content from the buffer. For compatibility with parent class
        :return:
        """
        if "w" in mode: self.reset_buffer(rel_path)
        buffer = self.save_buffer(rel_path)
        if not "a" in mode: buffer.seek(0)
        return buffer

    def reset_buffer(self, rel_path: str):
        rel_path = Path(self._format_relative_path(rel_path))
        self.file_dict.pop(rel_path, None)

    def register_file(self, rel_path: str, true_path:  Path = None):
        true_path = true_path or self.file_path(rel_path)
        rel_path = Path(self._format_relative_path(rel_path))
        return self.file_dict.setdefault(rel_path, true_path)

    def file_path(self, rel_path: str):
        if self.output_dir is None:
            raise UnsupportedOperation("Can't get a file path if output_dir is None")
        return super().file_path(rel_path)

    def to_zip(self, name) -> tuple[bytes, str]:
        archive_name = name.replace("[fname]", self.file_stem)
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for path, content in self.buffer_dict.items():
                zipf.writestr(zinfo_or_arcname=str(path), data=content.getvalue())
            for path, abs_path in self.file_dict.items():
                zipf.write(filename=abs_path, arcname=path)
        zip_buffer.seek(0)
        return zip_buffer.getvalue(), archive_name